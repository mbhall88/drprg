use rust_htslib::bcf;

const OGT_TAG: &str = "OGT";

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
    }
    // todo check for minor allele
    // todo rewrite gt
}

#[cfg(test)]
mod tests {
    use rust_htslib::bcf;
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
    }
}
