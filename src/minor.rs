use bstr::ByteSlice;
use drprg::VcfExt;
use rust_htslib::bcf;
use rust_htslib::bcf::record::GenotypeAllele;

const OGT_TAG: &str = "OGT";
const PDP_TAG: &str = "PDP";

pub struct MinorAllele {
    min_allele_freq: f32,
    max_gap: f32,
    max_gap_diff: f32,
}

impl MinorAllele {
    pub fn new(maf: f32, max_gap: f32, max_gap_diff: f32) -> Self {
        MinorAllele {
            min_allele_freq: maf,
            max_gap,
            max_gap_diff,
        }
    }
    pub fn add_vcf_headers(&self, header: &mut bcf::Header) {
        header.push_record(format!("##INFO=<ID={},Number=1,Type=String,Description=\"Original genotype after adjusting for minor allele depth proportions of {}\">", OGT_TAG, self.min_allele_freq).as_bytes());
        header.push_record(format!("##INFO=<ID={},Number=R,Type=Float,Description=\"Proportion of the total position depth found on this allele\">", PDP_TAG).as_bytes());
    }
    fn add_proportions_tag(
        record: &mut bcf::Record,
    ) -> Result<(), rust_htslib::errors::Error> {
        let pdp = record.depth_proportions();
        if let Some(v) = pdp {
            record.push_info_float(PDP_TAG.as_bytes(), &v)
        } else {
            Ok(())
        }
    }
    /// Check if any non-called allele has more than the min_allele_freq and return the index of the
    /// allele.
    /// Note: if multiple alleles over the min. then will return the one with the highest proportion
    /// Returns the index of the minor allele (if there is one) or 0 otherwise
    pub fn check_for_minor_alternate(
        &self,
        record: &mut bcf::Record,
    ) -> Result<isize, rust_htslib::errors::Error> {
        Self::add_proportions_tag(record)?;
        let dp_props = record.depth_proportions();

        let gt = record.called_allele();
        if record.allele_count() < 2 || dp_props.is_none() || gt < 0 {
            return Ok(-1);
        }

        let dp_props = dp_props.unwrap();
        let mut ix: Vec<(usize, &f32)> = dp_props.iter().enumerate().collect();
        ix.sort_by(|a, b| a.1.total_cmp(b.1));

        let mut largest_non_called: Option<(usize, f32)> = None;
        for (i, d) in ix.iter().rev() {
            if *i as i32 == gt {
                continue;
            } else {
                let gaps = record.format(b"GAPS").float()?[0][*i];
                let gaps_diff =
                    (gaps - record.format(b"GAPS").float()?[0][gt as usize]).abs();
                if **d >= self.min_allele_freq
                    && gaps <= self.max_gap
                    && gaps_diff <= self.max_gap_diff
                {
                    largest_non_called = Some((*i, **d));
                    break;
                }
            }
        }
        match largest_non_called {
            Some((i, _)) => Ok(i as isize),
            None => Ok(-1),
        }
    }
    pub fn adjust_genotype(
        record: &mut bcf::Record,
        new_gt: i32,
    ) -> Result<(), rust_htslib::errors::Error> {
        let og_gt = record.called_allele().to_string();
        record.push_info_string(OGT_TAG.as_bytes(), &[og_gt.as_bytes()])?;
        record.push_genotypes(&[GenotypeAllele::Unphased(new_gt)])
    }
    pub fn undo_genotype_adjustment(
        record: &mut bcf::Record,
    ) -> Result<(), rust_htslib::errors::Error> {
        let ogt = record.info(OGT_TAG.as_bytes()).string()?;
        if let Some(bb) = ogt {
            let original_gt_s = bb[0].to_str_lossy().to_string();
            let original_gt = original_gt_s.parse::<i32>().unwrap();
            record.push_genotypes(&[GenotypeAllele::Unphased(original_gt)])?;
            record.clear_info_string(OGT_TAG.as_bytes())?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use rust_htslib::bcf;
    use rust_htslib::bcf::record::GenotypeAllele;
    use tempfile::NamedTempFile;

    use super::*;

    #[test]
    fn test_write_vcf_headers() {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = bcf::Header::new();
        let ma = MinorAllele::new(1.0, 0.5, 0.1);
        ma.add_vcf_headers(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let view = vcf.header();

        assert!(view.name_to_id(OGT_TAG.as_bytes()).is_ok());
        assert!(view.name_to_id(PDP_TAG.as_bytes()).is_ok());
    }

    #[test]
    fn test_check_for_minor_alternate_null_call() {
        let ma = MinorAllele::new(0.5, 0.5, 0.0);
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = bcf::Header::new();

        header.push_sample(b"sample").push_record(br#"##FORMAT=<ID=MEAN_FWD_COVG,Number=R,Type=Integer,Description="Med forward coverage">"#).push_record(br#"##FORMAT=<ID=MEAN_REV_COVG,Number=R,Type=Integer,Description="Med reverse coverage">"#).push_record(br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#);
        ma.add_vcf_headers(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"A", b"T"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(-1)])
            .unwrap();
        record
            .push_format_integer(b"MEAN_FWD_COVG", &[5, 20])
            .expect("Failed to set forward coverage");
        record
            .push_format_integer(b"MEAN_REV_COVG", &[6, 30])
            .expect("Failed to set reverse coverage");

        let actual = ma.check_for_minor_alternate(&mut record).unwrap();
        let expected = -1;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_check_for_minor_alternate_alt_call() {
        let ma = MinorAllele::new(0.1, 0.5, 0.1);
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = bcf::Header::new();

        header.push_sample(b"sample").push_record(br#"##FORMAT=<ID=MEAN_FWD_COVG,Number=R,Type=Integer,Description="Med forward coverage">"#).push_record(br#"##FORMAT=<ID=MEAN_REV_COVG,Number=R,Type=Integer,Description="Med reverse coverage">"#).push_record(br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#).push_record(br#"##FORMAT=<ID=GAPS,Number=R,Type=Float,Description="Gaps">"#);
        ma.add_vcf_headers(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"A", b"T"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(1)])
            .unwrap();
        record
            .push_format_integer(b"MEAN_FWD_COVG", &[5, 20])
            .expect("Failed to set forward coverage");
        record
            .push_format_integer(b"MEAN_REV_COVG", &[6, 30])
            .expect("Failed to set reverse coverage");
        record
            .push_format_float(b"GAPS", &[0.0, 0.0])
            .expect("Failed to set gaps");

        let actual = ma.check_for_minor_alternate(&mut record).unwrap();
        let expected = 0;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_check_for_minor_alternate_ref_call_alt_has_most_depth() {
        let ma = MinorAllele::new(0.5, 0.5, 0.1);
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = bcf::Header::new();

        header.push_sample(b"sample").push_record(br#"##FORMAT=<ID=MEAN_FWD_COVG,Number=R,Type=Integer,Description="Med forward coverage">"#).push_record(br#"##FORMAT=<ID=MEAN_REV_COVG,Number=R,Type=Integer,Description="Med reverse coverage">"#).push_record(br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#).push_record(br#"##FORMAT=<ID=GAPS,Number=R,Type=Float,Description="Gaps">"#);
        ma.add_vcf_headers(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"A", b"T"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        record
            .push_format_integer(b"MEAN_FWD_COVG", &[5, 20])
            .expect("Failed to set forward coverage");
        record
            .push_format_integer(b"MEAN_REV_COVG", &[6, 30])
            .expect("Failed to set reverse coverage");
        record
            .push_format_float(b"GAPS", &[0.0, 0.0])
            .expect("Failed to set gaps");

        let actual = ma.check_for_minor_alternate(&mut record).unwrap();
        let expected = 1;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_check_for_minor_alternate_ref_call_ref_has_most_depth_alt_below_threshold()
    {
        let ma = MinorAllele::new(0.5, 0.5, 0.3);
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = bcf::Header::new();

        header.push_sample(b"sample").push_record(br#"##FORMAT=<ID=MEAN_FWD_COVG,Number=R,Type=Integer,Description="Med forward coverage">"#).push_record(br#"##FORMAT=<ID=MEAN_REV_COVG,Number=R,Type=Integer,Description="Med reverse coverage">"#).push_record(br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#).push_record(br#"##FORMAT=<ID=GAPS,Number=R,Type=Float,Description="Gaps">"#);
        ma.add_vcf_headers(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"A", b"T"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        record
            .push_format_integer(b"MEAN_FWD_COVG", &[50, 20])
            .expect("Failed to set forward coverage");
        record
            .push_format_integer(b"MEAN_REV_COVG", &[600, 30])
            .expect("Failed to set reverse coverage");
        record
            .push_format_float(b"GAPS", &[0.0, 0.2])
            .expect("Failed to set gaps");

        let actual = ma.check_for_minor_alternate(&mut record).unwrap();
        let expected = -1;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_check_for_minor_alternate_ref_call_ref_has_most_depth_alt_eq_threshold() {
        let ma = MinorAllele::new(50.0 / 160.0, 0.5, 0.1);
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = bcf::Header::new();

        header.push_sample(b"sample").push_record(br#"##FORMAT=<ID=MEAN_FWD_COVG,Number=R,Type=Integer,Description="Med forward coverage">"#).push_record(br#"##FORMAT=<ID=MEAN_REV_COVG,Number=R,Type=Integer,Description="Med reverse coverage">"#).push_record(br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#).push_record(br#"##FORMAT=<ID=GAPS,Number=R,Type=Float,Description="Gaps">"#);
        ma.add_vcf_headers(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"A", b"T"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        record
            .push_format_integer(b"MEAN_FWD_COVG", &[50, 20])
            .expect("Failed to set forward coverage");
        record
            .push_format_integer(b"MEAN_REV_COVG", &[60, 30])
            .expect("Failed to set reverse coverage");
        record
            .push_format_float(b"GAPS", &[0.0f32, 0.0f32])
            .expect("Failed to set gaps");

        let actual = ma.check_for_minor_alternate(&mut record).unwrap();
        let expected = 1;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_check_for_minor_alternate_ref_call_ref_has_most_depth_alt_above_threshold()
    {
        let ma = MinorAllele::new(50.0 / 160.0, 0.5, 0.1);
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = bcf::Header::new();

        header.push_sample(b"sample").push_record(br#"##FORMAT=<ID=MEAN_FWD_COVG,Number=R,Type=Integer,Description="Med forward coverage">"#).push_record(br#"##FORMAT=<ID=MEAN_REV_COVG,Number=R,Type=Integer,Description="Med reverse coverage">"#).push_record(br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#).push_record(br#"##FORMAT=<ID=GAPS,Number=R,Type=Float,Description="Gaps">"#);
        ma.add_vcf_headers(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"A", b"T"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        record
            .push_format_integer(b"MEAN_FWD_COVG", &[50, 21])
            .expect("Failed to set forward coverage");
        record
            .push_format_integer(b"MEAN_REV_COVG", &[60, 30])
            .expect("Failed to set reverse coverage");
        record
            .push_format_float(b"GAPS", &[0.0, 0.0])
            .expect("Failed to set gaps");

        let actual = ma.check_for_minor_alternate(&mut record).unwrap();
        let expected = 1;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_check_for_minor_alternate_ref_call_ref_has_most_depth_alt_below_gaps_threshold(
    ) {
        let ma = MinorAllele::new(50.0 / 160.0, 0.4, 0.5);
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = bcf::Header::new();

        header.push_sample(b"sample").push_record(br#"##FORMAT=<ID=MEAN_FWD_COVG,Number=R,Type=Integer,Description="Med forward coverage">"#).push_record(br#"##FORMAT=<ID=MEAN_REV_COVG,Number=R,Type=Integer,Description="Med reverse coverage">"#).push_record(br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#).push_record(br#"##FORMAT=<ID=GAPS,Number=R,Type=Float,Description="Gaps">"#);
        ma.add_vcf_headers(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"A", b"T"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        record
            .push_format_integer(b"MEAN_FWD_COVG", &[50, 21])
            .expect("Failed to set forward coverage");
        record
            .push_format_integer(b"MEAN_REV_COVG", &[60, 30])
            .expect("Failed to set reverse coverage");
        record
            .push_format_float(b"GAPS", &[0.0, 0.45])
            .expect("Failed to set gaps");

        let actual = ma.check_for_minor_alternate(&mut record).unwrap();
        let expected = -1;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_check_for_minor_alternate_ref_call_no_depth() {
        let ma = MinorAllele::new(0.1, 0.5, 0.0);
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = bcf::Header::new();

        header.push_sample(b"sample").push_record(br#"##FORMAT=<ID=MEAN_FWD_COVG,Number=R,Type=Integer,Description="Med forward coverage">"#).push_record(br#"##FORMAT=<ID=MEAN_REV_COVG,Number=R,Type=Integer,Description="Med reverse coverage">"#).push_record(br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#);
        ma.add_vcf_headers(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"A", b"T"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        record
            .push_format_integer(b"MEAN_FWD_COVG", &[0, 0])
            .expect("Failed to set forward coverage");
        record
            .push_format_integer(b"MEAN_REV_COVG", &[0, 0])
            .expect("Failed to set reverse coverage");

        let actual = ma.check_for_minor_alternate(&mut record).unwrap();
        let expected = -1;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_check_for_minor_alternate_calls_alternate_but_other_alt_is_minor() {
        let ma = MinorAllele::new(0.2, 0.3, 0.1);
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = bcf::Header::new();

        header.push_sample(b"sample").push_record(br#"##FORMAT=<ID=MEAN_FWD_COVG,Number=R,Type=Integer,Description="Med forward coverage">"#).push_record(br#"##FORMAT=<ID=MEAN_REV_COVG,Number=R,Type=Integer,Description="Med reverse coverage">"#).push_record(br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#).push_record(br#"##FORMAT=<ID=GAPS,Number=R,Type=Float,Description="Gaps">"#);
        ma.add_vcf_headers(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"A", b"T", b"C", b"G"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(3)])
            .unwrap();
        record
            .push_format_integer(b"MEAN_FWD_COVG", &[0, 21, 1, 70])
            .expect("Failed to set forward coverage");
        record
            .push_format_integer(b"MEAN_REV_COVG", &[1, 30, 0, 70])
            .expect("Failed to set reverse coverage");
        record
            .push_format_float(b"GAPS", &[1.0, 0.0, 1.0, 0.0])
            .expect("Failed to set gaps");

        let actual = ma.check_for_minor_alternate(&mut record).unwrap();
        let expected = 1;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_check_for_minor_alternate_below_threshold_but_above_diff() {
        let ma = MinorAllele::new(50.0 / 160.0, 0.4, 0.1);
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = bcf::Header::new();

        header.push_sample(b"sample").push_record(br#"##FORMAT=<ID=MEAN_FWD_COVG,Number=R,Type=Integer,Description="Med forward coverage">"#).push_record(br#"##FORMAT=<ID=MEAN_REV_COVG,Number=R,Type=Integer,Description="Med reverse coverage">"#).push_record(br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#).push_record(br#"##FORMAT=<ID=GAPS,Number=R,Type=Float,Description="Gaps">"#);
        ma.add_vcf_headers(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"A", b"T"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        record
            .push_format_integer(b"MEAN_FWD_COVG", &[50, 21])
            .expect("Failed to set forward coverage");
        record
            .push_format_integer(b"MEAN_REV_COVG", &[60, 30])
            .expect("Failed to set reverse coverage");
        record
            .push_format_float(b"GAPS", &[0.0f32, 0.25f32])
            .expect("Failed to set gaps");

        let actual = ma.check_for_minor_alternate(&mut record).unwrap();
        let expected = -1;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_check_for_minor_alternate_above_threshold_below_diff() {
        let ma = MinorAllele::new(50.0 / 160.0, 0.4, 0.1);
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();
        let mut header = bcf::Header::new();

        header.push_sample(b"sample").push_record(br#"##FORMAT=<ID=MEAN_FWD_COVG,Number=R,Type=Integer,Description="Med forward coverage">"#).push_record(br#"##FORMAT=<ID=MEAN_REV_COVG,Number=R,Type=Integer,Description="Med reverse coverage">"#).push_record(br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#).push_record(br#"##FORMAT=<ID=GAPS,Number=R,Type=Float,Description="Gaps">"#);
        ma.add_vcf_headers(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = vcf.empty_record();
        let alleles: &[&[u8]] = &[b"A", b"T"];
        record.set_alleles(alleles).expect("Failed to set alleles");
        record
            .push_genotypes(&[GenotypeAllele::Unphased(0)])
            .unwrap();
        record
            .push_format_integer(b"MEAN_FWD_COVG", &[50, 21])
            .expect("Failed to set forward coverage");
        record
            .push_format_integer(b"MEAN_REV_COVG", &[60, 30])
            .expect("Failed to set reverse coverage");
        record
            .push_format_float(b"GAPS", &[0.39f32, 0.45f32])
            .expect("Failed to set gaps");

        let actual = ma.check_for_minor_alternate(&mut record).unwrap();
        let expected = -1;

        assert_eq!(actual, expected)
    }
}
