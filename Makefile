# path macros
EXTDIR := src/ext

# URLs
PANDORA_URL := "https://github.com/rmcolq/pandora/releases/download/0.9.0-rc2/pandora-linux-precompiled-v0.9.0-rc2"
MAKEPRG_URL := "https://github.com/leoisl/make_prg/releases/download/v0.2.0_prototype/make_prg_0.2.0_prototype"

# binary names
PANDORA := $(EXTDIR)/pandora
MAKEPRG := $(EXTDIR)/make_prg

# clean files list
CLEAN_LIST := $(PANDORA) $(MAKEPRG)

define download
      wget $(1) -O $(2)
      chmod +x $(2)
endef

.PHONY: makedir
makedir:
	@mkdir -p $(EXTDIR)

.PHONY: pandora
pandora: makedir
	$(call download,$(PANDORA_URL),$(PANDORA))

.PHONY: makeprg
makeprg: makedir
	$(call download,$(MAKEPRG_URL),$(MAKEPRG))

.PHONY: clean
clean:
	@echo CLEAN $(CLEAN_LIST)
	@rm -f $(CLEAN_LIST)

.PHONY: deps
deps: pandora makeprg