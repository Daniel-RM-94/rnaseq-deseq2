#!/usr/bin/env bash
# Comprueba la contraseña de RStudio Server antes de arrancarlo.
# Si no hay una contraseña segura, el contenedor no se inicia.
set -euo pipefail

readonly MIN_LENGTH=12
readonly WEAK_PASSWORDS=(
  "rstudio" "bioconductor" "password" "changeme" "cambia-esta-contraseña"
)

fail() {
  printf '\nERROR: %s\n' "$1" >&2
  printf 'Set RSTUDIO_PASSWORD in the .env file (see .env.example).\n\n' >&2
  exit 1
}

if [[ "${DISABLE_AUTH:-}" == "true" ]]; then
  fail "DISABLE_AUTH is not allowed: RStudio Server must require authentication."
fi

password="${PASSWORD:-}"

if [[ -z "$password" ]]; then
  fail "No RStudio Server password was provided."
fi

if (( ${#password} < MIN_LENGTH )); then
  fail "The RStudio Server password must be at least ${MIN_LENGTH} characters long."
fi

for weak in "${WEAK_PASSWORDS[@]}"; do
  if [[ "${password,,}" == "$weak" ]]; then
    fail "The RStudio Server password is a default or example value."
  fi
done

unset password
exec /init "$@"
