# path macros
EXTDIR := $(abspath src/ext)
MAFFT_SRC := $(EXTDIR)/mafft

# URLs
PANDORA_URL := "https://github.com/rmcolq/pandora/releases/download/0.9.0/pandora-linux-precompiled-v0.9.0"
MAKEPRG_URL := "https://github.com/leoisl/make_prg/releases/download/v0.2.0/make_prg_0.2.0"
BCFTOOLS_URL := "https://github.com/samtools/bcftools/releases/download/1.12/bcftools-1.12.tar.bz2"
MAFFT_URL := "https://mafft.cbrc.jp/alignment/software/mafft-7.475-without-extensions-src.tgz"

# binary names
PANDORA := $(EXTDIR)/pandora
MAKEPRG := $(EXTDIR)/make_prg
MAFFT := $(EXTDIR)/mafft/bin/mafft
BCFTOOLS := $(EXTDIR)/bcftools

# clean files list
CLEAN_LIST := $(PANDORA) $(MAKEPRG) $(MAFFT)

define download
      wget $(1) -O $(2)
      chmod +x $(2)
endef

.PHONY: makedir
makedir:
	@mkdir -p $(EXTDIR)

.PHONY: bcftools
bcftools: makedir
	cd $(EXTDIR) \
	&& ( wget $(BCFTOOLS_URL) -O - | tar -xjf - ) \
	&& cd bcftools-1.12 \
	&& ./configure --prefix=$(EXTDIR)/bcftools-1.12 \
	&& make \
	&& make install \
	&& cd .. \
	&& ln -s bcftools-1.12/bin/bcftools $(BCFTOOLS)

.PHONY: mafft
mafft: makedir
	mkdir -p $(MAFFT_SRC)
	wget -O - $(MAFFT_URL) | tar xzf - -C $(MAFFT_SRC) --strip-components=1
	sed -i '1s?/usr/local?$(MAFFT_SRC)?' $(MAFFT_SRC)/core/Makefile
	make -C $(MAFFT_SRC)/core install PREFIX=$(MAFFT_SRC)

.PHONY: pandora
pandora: makedir
	$(call download,$(PANDORA_URL),$(PANDORA))

.PHONY: makeprg
makeprg: makedir
	$(call download,$(MAKEPRG_URL),$(MAKEPRG))

.PHONY: clean
clean:
	@echo CLEAN $(CLEAN_LIST)
	@rm -rf $(CLEAN_LIST)

.PHONY: deps
deps: pandora makeprg mafft bcftools