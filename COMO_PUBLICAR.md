# Cómo subir este taller a GitHub Pages

## 1. Crear el repositorio

1. Entra a https://github.com/new e inicia sesión.
2. Escribe `Calculo-Vectorial-Taller-1` como nombre (o utiliza otro nombre).
3. Selecciona **Public** para usar GitHub Pages con GitHub Free.
4. Crea el repositorio. Puedes dejar sin marcar la creación de README, porque el ZIP ya incluye uno.

## 2. Subir los archivos

1. **Descomprime** el ZIP que acompaña a este proyecto.
2. En la página del repositorio, elige **Add file → Upload files**.
3. Arrastra o selecciona **todos los archivos de dentro de la carpeta** `Calculo_Vectorial_Taller_1_Pluto_Oscuro` (index.html, los 18 archivos Parte*.html/jl, README.md, COMO_PUBLICAR.md, favicon.svg, pluto-site-theme.css, pluto-site-theme.js y .nojekyll). En algunos sistemas los archivos que empiezan con punto están ocultos; si no puedes seleccionar `.nojekyll`, los demás archivos también pueden funcionar, pero conviene incluirlo.
4. Confirma la subida con **Commit changes**.

**Importante:** No subas únicamente `t1.rar` ni el archivo ZIP; GitHub Pages necesita los HTML descomprimidos. `index.html` debe estar en la raíz, al mismo nivel que `Parte1.html`.

### Alternativa: con Git instalado

Desde dentro de la carpeta del proyecto:

```bash
git init
git add .
git commit -m "Publicar Taller 1 de Cálculo Vectorial"
git branch -M main
git remote add origin https://github.com/TUUSUARIO/Calculo-Vectorial-Taller-1.git
git push -u origin main
```

Crea primero el repositorio remoto vacío y reemplaza `TUUSUARIO`.

## 3. Activar GitHub Pages

1. Abre el repositorio en GitHub y entra en **Settings**.
2. En el menú lateral entra en **Pages**.
3. En **Build and deployment → Source**, selecciona **Deploy from a branch**.
4. Selecciona la rama **main** y la carpeta **/(root)**.
5. Pulsa **Save** y comprueba que aparezca la URL de tu sitio.

### Direcciones de ejemplo

- Inicio: `https://TUUSUARIO.github.io/Calculo-Vectorial-Taller-1/`
- Parte 1: `https://TUUSUARIO.github.io/Calculo-Vectorial-Taller-1/Parte1.html`
- Parte 9: `https://TUUSUARIO.github.io/Calculo-Vectorial-Taller-1/Parte9.html`

Si tu repositorio tiene otro nombre, cámbialo en esas direcciones.

## 4. Comprobar que quedó bien

- Abre la dirección principal y verifica que aparecen las nueve partes.
- Usa el campo de búsqueda y cambia entre modo claro y oscuro.
- Abre la Parte 1 y la Parte 9 para confirmar que cargan las fórmulas y resultados y que el modo oscuro se mantiene al navegar. El icono de sol permite cambiar a modo claro.
- Comprueba que el enlace **Ver código .jl** abre el archivo original.

## Notas técnicas

- El índice es HTML, CSS y JavaScript puro. No requiere React, npm, Hugo ni compilación.
- Los nueve HTML conservan el contenido original de Pluto, incluidos los recursos embebidos. Se añadió una barra flotante de retorno al índice y dos archivos compartidos de tema visual.
- Las exportaciones HTML permiten **leer resultados guardados**, pero no ejecutar de nuevo los cálculos Julia dentro de GitHub Pages.
- Los nueve notebooks originales `.jl` se incluyen para estudio y reutilización.
- El repositorio y la web resultante serán públicos: revisa el contenido antes de hacer la subida.

**Referencias oficiales:**
- GitHub Pages: https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site
- Exportación de Pluto: https://plutojl.org/en/docs/export-html/
