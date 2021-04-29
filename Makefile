# path macros
EXTDIR := $(abspath src/ext)
MAFFT_SRC := $(EXTDIR)/mafft

# URLs
PANDORA_URL := "https://github.com/rmcolq/pandora/releases/download/0.9.0-rc2/pandora-linux-precompiled-v0.9.0-rc2"
MAKEPRG_URL := "https://github.com/leoisl/make_prg/releases/download/v0.2.0_prototype/make_prg_0.2.0_prototype"
MAFFT_URL := "https://mafft.cbrc.jp/alignment/software/mafft-7.475-without-extensions-src.tgz"

# binary names
PANDORA := $(EXTDIR)/pandora
MAKEPRG := $(EXTDIR)/make_prg
MAFFT := $(EXTDIR)/mafft/bin/mafft

# clean files list
CLEAN_LIST := $(PANDORA) $(MAKEPRG) $(MAFFT)

define download
      wget $(1) -O $(2)
      chmod +x $(2)
endef

.PHONY: makedir
makedir:
	@mkdir -p $(EXTDIR)

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
deps: pandora makeprg mafft