use bstr::ByteSlice;
use clap::builder::ArgPredicate;
use clap::Parser;
use drprg::VcfExt;
use float_cmp::approx_eq;
use rust_htslib::bcf;
use rust_htslib::bcf::record::GenotypeAllele;

const OGT_TAG: &str = "OGT";
const PDP_TAG: &str = "PDP";
pub static MINOR_AF: f32 = 1.0;
pub static MINOR_AF_ILLUMINA: &str = "0.1";
pub static MAX_GAPS: f32 = 0.5;
pub static MAX_CALLED_GAPS: f32 = 0.39;
pub static MAX_GAPS_DIFF: f32 = 0.2;
pub const MINOR_MIN_COVG: i32 = 3;
pub const MINOR_MIN_STRAND_BIAS: f32 = 0.01;

#[derive(Parser, Debug, Default)]
pub struct MinorAllele {
    /// Minimum allele frequency to call variants
    ///
    /// If an alternate allele has at least this fraction of the depth, a minor resistance
    /// ("r") prediction is made. Set to 1 to disable. If --illumina is passed, the default
    /// is 0.1
    #[clap(
    short = 'f',
    long,
    value_name = "FLOAT[0.0-1.0]",
    default_value_t = MINOR_AF,
    default_value_if("is_illumina", ArgPredicate::IsPresent, MINOR_AF_ILLUMINA)
    )]
    pub(crate) maf: f32,
    /// Maximum value for the GAPS tag when calling a minor allele
    #[clap(long, default_value_t = MAX_GAPS, hide = true)]
    pub(crate) max_gaps: f32,
    /// Maximum value for the GAPS tag of the called allele when calling a minor allele
    #[clap(long, default_value_t = MAX_CALLED_GAPS, hide = true)]
    pub(crate) max_called_gaps: f32,
    /// Maximum difference allowed between major and minor alleles for the GAPS tag when calling a minor allele
    #[clap(long, default_value_t = MAX_GAPS_DIFF, hide = true)]
    pub(crate) max_gaps_diff: f32,
    /// Minimum depth allowed on a minor allele
    #[clap(long, default_value_t = MINOR_MIN_COVG, hide = true)]
    pub(crate) minor_min_covg: i32,
    /// Minimum strand bias ratio allowed on a minor allele
    #[clap(long, default_value_t = MINOR_MIN_STRAND_BIAS, hide = true)]
    pub(crate) minor_min_strand_bias: f32,
}

impl MinorAllele {
    pub fn add_vcf_headers(&self, header: &mut bcf::Header) {
        header.push_record(format!("##INFO=<ID={},Number=1,Type=String,Description=\"Original genotype after adjusting for minor allele depth proportions of {}\">", OGT_TAG, self.maf).as_bytes());
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
        let called_gaps = record.format(b"GAPS").float()?[0][gt as usize];

        if called_gaps > self.max_called_gaps {
            return Ok(-1);
        }

        let mut largest_non_called: Option<(usize, f32)> = None;
        for (i, d) in ix.iter().rev() {
            if *i as i32 == gt {
                continue;
            } else {
                let gaps = record.format(b"GAPS").float()?[0][*i];
                let gaps_diff = gaps - record.format(b"GAPS").float()?[0][gt as usize];
                if **d >= self.maf
                    && gaps <= self.max_gaps
                    && gaps_diff <= self.max_gaps_diff
                {
                    largest_non_called = Some((*i, **d));
                    break;
                }
            }
        }
        match largest_non_called {
            Some((gt, _)) => {
                let (fc, rc) = record.coverage().unwrap_or_else(|| (vec![0], vec![0]));
                let sum_covg = (fc[gt] + rc[gt]) as f32;
                let covg = fc.get(gt).unwrap_or(&0) + rc.get(gt).unwrap_or(&0);
                let has_low_covg = covg < self.minor_min_covg;
                let has_strand_bias = if approx_eq!(f32, sum_covg, 0.0) {
                    true
                } else {
                    fc[gt].min(rc[gt]) as f32 / sum_covg < self.minor_min_strand_bias
                };

                if has_low_covg || has_strand_bias {
                    Ok(-1)
                } else {
                    Ok(gt as isize)
                }
            }
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
        let ma = MinorAllele {
            maf: 1.0,
            max_gaps: 0.5,
            max_called_gaps: 0.5,
            max_gaps_diff: 0.1,
            minor_min_covg: 0,
            minor_min_strand_bias: 0.0,
        };
        ma.add_vcf_headers(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let view = vcf.header();

        assert!(view.name_to_id(OGT_TAG.as_bytes()).is_ok());
        assert!(view.name_to_id(PDP_TAG.as_bytes()).is_ok());
    }

    #[test]
    fn test_check_for_minor_alternate_null_call() {
        let ma = MinorAllele {
            maf: 0.5,
            max_gaps: 0.5,
            max_called_gaps: 0.5,
            ..Default::default()
        };
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
        let ma = MinorAllele {
            maf: 0.1,
            max_gaps: 0.5,
            max_called_gaps: 0.5,
            max_gaps_diff: 0.1,
            ..Default::default()
        };
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
        let ma = MinorAllele {
            maf: 0.5,
            max_gaps: 0.5,
            max_called_gaps: 0.5,
            max_gaps_diff: 0.1,
            ..Default::default()
        };
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
        let ma = MinorAllele {
            maf: 0.5,
            max_gaps: 0.5,
            max_called_gaps: 0.5,
            max_gaps_diff: 0.3,
            minor_min_covg: 0,
            minor_min_strand_bias: 0.0,
        };
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
        let ma = MinorAllele {
            maf: 50.0 / 160.0,
            max_gaps: 0.5,
            max_called_gaps: 0.5,
            max_gaps_diff: 0.1,
            minor_min_covg: 0,
            minor_min_strand_bias: 0.0,
        };
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
        let ma = MinorAllele {
            maf: 50.0 / 160.0,
            max_gaps: 0.5,
            max_called_gaps: 0.5,
            max_gaps_diff: 0.1,
            minor_min_covg: 0,
            minor_min_strand_bias: 0.0,
        };
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
        let ma = MinorAllele {
            maf: 50.0 / 160.0,
            max_gaps: 0.4,
            max_called_gaps: 0.4,
            max_gaps_diff: 0.5,
            minor_min_covg: 0,
            minor_min_strand_bias: 0.0,
        };
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
        let ma = MinorAllele {
            maf: 0.1,
            max_gaps: 0.5,
            max_called_gaps: 0.5,
            max_gaps_diff: 0.0,
            minor_min_covg: 0,
            minor_min_strand_bias: 0.0,
        };
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
        let ma = MinorAllele {
            maf: 0.2,
            max_gaps: 0.3,
            max_called_gaps: 0.3,
            max_gaps_diff: 0.1,
            minor_min_covg: 0,
            minor_min_strand_bias: 0.0,
        };
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
        let ma = MinorAllele {
            maf: 50.0 / 160.0,
            max_gaps: 0.4,
            max_called_gaps: 0.4,
            max_gaps_diff: 0.1,
            minor_min_covg: 0,
            minor_min_strand_bias: 0.0,
        };
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
        let ma = MinorAllele {
            maf: 50.0 / 160.0,
            max_gaps: 0.4,
            max_called_gaps: 0.4,
            max_gaps_diff: 0.1,
            minor_min_covg: 0,
            minor_min_strand_bias: 0.0,
        };
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

    #[test]
    fn test_check_for_minor_alternate_alt_has_less_gaps_than_ref() {
        let ma = MinorAllele {
            maf: 0.1,
            max_gaps: 0.4,
            max_called_gaps: 0.4,
            max_gaps_diff: 0.1,
            minor_min_covg: 0,
            minor_min_strand_bias: 0.0,
        };
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
            .push_format_integer(b"MEAN_FWD_COVG", &[64, 13])
            .expect("Failed to set forward coverage");
        record
            .push_format_integer(b"MEAN_REV_COVG", &[50, 12])
            .expect("Failed to set reverse coverage");
        record
            .push_format_float(b"GAPS", &[0.3333, 0.0])
            .expect("Failed to set gaps");

        let actual = ma.check_for_minor_alternate(&mut record).unwrap();
        let expected = 1;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_check_for_minor_alternate_low_covg() {
        let ma = MinorAllele {
            maf: 0.1,
            max_gaps: 0.3,
            max_called_gaps: 0.3,
            max_gaps_diff: 0.1,
            minor_min_covg: 3,
            minor_min_strand_bias: 0.0,
        };
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
            .push_format_integer(b"MEAN_FWD_COVG", &[6, 1])
            .expect("Failed to set forward coverage");
        record
            .push_format_integer(b"MEAN_REV_COVG", &[5, 1])
            .expect("Failed to set reverse coverage");
        record
            .push_format_float(b"GAPS", &[0.3333, 0.0])
            .expect("Failed to set gaps");

        let actual = ma.check_for_minor_alternate(&mut record).unwrap();
        let expected = -1;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_check_for_minor_alternate_low_stand_bias() {
        let ma = MinorAllele {
            maf: 0.1,
            max_gaps: 0.3,
            max_called_gaps: 0.3,
            max_gaps_diff: 0.1,
            minor_min_covg: 3,
            minor_min_strand_bias: 0.01,
        };
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
            .push_format_integer(b"MEAN_FWD_COVG", &[6, 3])
            .expect("Failed to set forward coverage");
        record
            .push_format_integer(b"MEAN_REV_COVG", &[5, 0])
            .expect("Failed to set reverse coverage");
        record
            .push_format_float(b"GAPS", &[0.3333, 0.0])
            .expect("Failed to set gaps");

        let actual = ma.check_for_minor_alternate(&mut record).unwrap();
        let expected = -1;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_check_for_minor_alternate_low_stand_bias_and_covg() {
        let ma = MinorAllele {
            maf: 0.1,
            max_gaps: 0.3,
            max_called_gaps: 0.3,
            max_gaps_diff: 0.1,
            minor_min_covg: 3,
            minor_min_strand_bias: 0.01,
        };
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
            .push_format_integer(b"MEAN_FWD_COVG", &[6, 2])
            .expect("Failed to set forward coverage");
        record
            .push_format_integer(b"MEAN_REV_COVG", &[5, 0])
            .expect("Failed to set reverse coverage");
        record
            .push_format_float(b"GAPS", &[0.3333, 0.0])
            .expect("Failed to set gaps");

        let actual = ma.check_for_minor_alternate(&mut record).unwrap();
        let expected = -1;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_check_for_minor_alternate_called_allele_over_max_called_gap() {
        let ma = MinorAllele {
            maf: 0.1,
            max_gaps: 0.5,
            max_called_gaps: 0.39,
            max_gaps_diff: 0.2,
            minor_min_covg: 3,
            minor_min_strand_bias: 0.01,
        };
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
            .push_format_integer(b"MEAN_FWD_COVG", &[20, 16])
            .expect("Failed to set forward coverage");
        record
            .push_format_integer(b"MEAN_REV_COVG", &[11, 8])
            .expect("Failed to set reverse coverage");
        record
            .push_format_float(b"GAPS", &[0.4, 0.5])
            .expect("Failed to set gaps");

        let actual = ma.check_for_minor_alternate(&mut record).unwrap();
        let expected = -1;

        assert_eq!(actual, expected)
    }
}
