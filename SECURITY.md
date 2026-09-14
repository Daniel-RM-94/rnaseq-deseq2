# Política de seguridad

## Reportar una vulnerabilidad

**No abras un issue público** para reportar un problema de seguridad. Usa la
opción *Report a vulnerability* de la pestaña **Security** del repositorio
(GitHub Private Vulnerability Reporting). Los reportes se atienden de forma
privada y se publican solo cuando exista una solución.

## Gestión de credenciales

- **Sin credenciales en el código.** El repositorio no contiene contraseñas,
  tokens ni valores por defecto.
- **RStudio Server exige una contraseña propia**, definida en `RSTUDIO_PASSWORD`
  dentro de `.env`. Ese archivo:
  - está excluido de git (`.gitignore`);
  - está excluido de la imagen Docker (`.dockerignore`).
- **Validación al arrancar** (`docker/start_rstudio.sh`). El contenedor no se
  inicia si:
  - la contraseña está vacía;
  - tiene menos de 12 caracteres;
  - coincide con un valor por defecto o de ejemplo;
  - se intenta desactivar la autenticación (`DISABLE_AUTH=true`).
- **Acceso solo local.** El puerto 8787 está enlazado a `127.0.0.1`. No lo
  cambies a `0.0.0.0`. Para acceder desde otro equipo, usa un túnel SSH:
  `ssh -L 8787:localhost:8787 usuario@servidor`.
- **CI sin secretos.** El pipeline y los tests no necesitan credenciales. Si en
  el futuro se necesitan, deben guardarse en *GitHub Actions Secrets*, nunca
  en el código ni en los workflows.

## Configuración recomendada del repositorio en GitHub

En *Settings → Code security*:

- [ ] **Secret scanning** y **Push protection**: bloquean commits que contengan
      secretos.
- [ ] **Dependabot alerts**: avisos de vulnerabilidades. `dependabot.yml`, ya
      incluido, mantiene actualizadas las GitHub Actions y la imagen base.
- [ ] **Private vulnerability reporting**.

En *Settings → Branches*:

- [ ] **Protección de `main`**: exigir *pull request* y que la CI pase antes de
      fusionar.

## Si se publica una credencial por error

1. **Considérala comprometida y cámbiala de inmediato.** Eliminarla del código
   no basta, porque sigue en el historial de git y puede haber sido copiada.
2. Elimínala del historial con
   [`git filter-repo`](https://github.com/newren/git-filter-repo) y fuerza la
   actualización del repositorio remoto.
3. Revisa los accesos que pudieran haberse producido con esa credencial.
