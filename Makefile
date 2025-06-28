PLINK2_LINUX_URL = https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20230726.zip
PLINK2_DARWIN_URL = https://s3.amazonaws.com/plink2-assets/plink2_mac_arm64_20240806.zip

PLINK2 = bin/plink2

UNAME_S := $(shell uname -s)
ARCH := $(shell uname -m)

.PHONY: all clean install_plink1 install_plink2 install 
.EXPORT_ALL_VARIABLES:

all:
	$(MAKE) install
	$(MAKE) install_plink2

activate:
	source ./.venv/bin/activate

install-dev:
	@echo "Installing development dependencies..."
	uv sync --extra dev

install-test:
	@echo "Installing test dependencies..."
	uv sync --extra test

install-all:
	@echo "Installing all dependencies..."
	uv sync --extra test --extra dev
	uv run pre-commit install

install:
	make install-all

install_plink2:
	@if [ "$(UNAME_S)" = "Linux" ]; then \
		echo "Detected Linux $(ARCH)"; \
		echo "Downloading and installing PLINK 2.0..."; \
		curl -LO $(PLINK2_LINUX_URL); \
		unzip plink2_linux_x86_64_20230726.zip; \
		mv plink2 $(PLINK2); \
		chmod +x $(PLINK2); \
		rm plink2_linux_x86_64_20230726.zip; \
		rm -f prettify; \
		rm -f toy.map; \
		rm -f toy.ped; \
		rm -f LICENSE; \
		echo "PLINK 2.0 installation completed."; \
	elif [ "$(UNAME_S)" = "Darwin" ]; then \
		echo "Detected macOS $(ARCH)"; \
		echo "Downloading and installing PLINK 2.0..."; \
		curl -LO $(PLINK2_DARWIN_URL); \
		unzip plink2_mac_arm64_20240806.zip; \
		mv plink2 $(PLINK2); \
		chmod +x $(PLINK2); \
		rm plink2_mac_arm64_20240806.zip; \
		rm -f prettify; \
		rm -f toy.map; \
		rm -f toy.ped; \
		rm -f LICENSE; \
		echo "PLINK 2.0 installation completed."; \
	else \
		echo "Unsupported platform: $(UNAME_S) $(ARCH)"; \
		exit 1; \
	fi

clean:
	@echo "Cleaning up..."
	@rm -f $(PLINK2)
	@echo "Clean up completed."

build-image:
	docker pull ubuntu