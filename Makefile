COMPOSE := docker compose

.DEFAULT_GOAL := help
.PHONY: help build pipeline test rstudio stop shell clean

help: ## Muestra esta ayuda
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-10s\033[0m %s\n", $$1, $$2}'

build: ## Construye la imagen Docker
	$(COMPOSE) build rstudio

pipeline: ## Ejecuta el análisis completo y genera el informe HTML
	$(COMPOSE) run --rm pipeline

test: ## Ejecuta los tests unitarios
	$(COMPOSE) run --rm tests

rstudio: ## Inicia RStudio Server (puerto RSTUDIO_PORT de .env, 8787 por defecto)
	@test -f .env || { echo "Falta .env: ejecuta 'cp .env.example .env' y define RSTUDIO_PASSWORD"; exit 1; }
	$(COMPOSE) up -d --wait rstudio
	@echo "RStudio Server disponible en http://$$($(COMPOSE) port rstudio 8787) (usuario: rstudio)"

stop: ## Detiene los contenedores
	$(COMPOSE) down

shell: ## Abre una terminal dentro del contenedor
	$(COMPOSE) run --rm pipeline bash

clean: ## Elimina los resultados y datos intermedios generados
	$(COMPOSE) run --rm pipeline find results data/processed data/raw -type f ! -name README.md ! -name .gitkeep -delete
