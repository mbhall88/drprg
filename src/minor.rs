use rust_htslib::bcf;
use drprg::VcfExt;

const OGT_TAG: &str = "OGT";
const PDP_TAG: &str = "PDP";

pub struct MinorAllele {
    min_allele_freq: f32,
}

impl MinorAllele {
    pub fn new(maf: f32) -> Self {
        MinorAllele {
            min_allele_freq: maf,
        }
    }
    pub fn add_vcf_headers(&self, header: &mut bcf::Header) {
        header.push_record(format!("##INFO=<ID={},Number=1,Type=String,Description=\"Original genotype after adjusting for minor allele depth proportions of {}\">", OGT_TAG, self.min_allele_freq).as_bytes());
        header.push_record(format!("##INFO=<ID={},Number=R,Type=Float,Description=\"Proportion of the total position depth found on this allele\">", PDP_TAG).as_bytes());
    }
    fn add_proportions_tag(record: &mut bcf::Record) -> Result<(), rust_htslib::errors::Error>{
        let pdp = record.depth_proportions();
        if let Some(v) = pdp {
            record.push_info_float(PDP_TAG.as_bytes(), &v)
        } else {
            Ok(())
        }
    }
    /// Check if any non-called allele has more than the min_allele_freq and return the index of the
    /// allele. Only applicable if record is called ref.
    /// Note: if multiple alleles over the min. then will return the one with the highest proportion
    pub fn check_for_minor_alternate(&self, record: &mut bcf::Record) -> Result<usize, rust_htslib::errors::Error> {
        Self::add_proportions_tag(record)?;
        let gt = record.called_allele();
        if gt != 0 || record.allele_count() < 2 {
            return Ok(0)
        }
        todo!()
    }
    // todo rewrite gt
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
        let ma = MinorAllele::new(1.0);
        ma.add_vcf_headers(&mut header);
        let vcf =
            bcf::Writer::from_path(path, &header, true, bcf::Format::Vcf).unwrap();
        let view = vcf.header();

        assert!(view.name_to_id(OGT_TAG.as_bytes()).is_ok());
        assert!(view.name_to_id(PDP_TAG.as_bytes()).is_ok());
    }
}
