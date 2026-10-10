# Taller 1 — Cálculo Vectorial

### [VER PÁGINA WEB — Abrir los 9 notebooks de Pluto.jl](https://jostinbb.github.io/Calculo-Vectorial/)

**Sitio web:** https://jostinbb.github.io/Calculo-Vectorial/  
**Autor:** Jostin Balarezo  

---

**Autor:** Jostin Balarezo  
**Carrera:** Ingeniería en Ciencias de la Computación — Universidad de Cuenca  
**Herramientas:** Julia, Pluto.jl, SymPy y otras bibliotecas usadas por los notebooks.

Este repositorio reúne nueve partes del Taller 1. La página principal `index.html`, diseñada con la estética de Pluto.jl y modo oscuro predeterminado, permite consultar cada exportación HTML y acceder al código original en Julia.

## Contenido

| Parte | Tema del notebook | HTML | Fuente Julia |
|:--:|---|---|---|
| 1 | Integrales de línea: función escalar | [Abrir](Parte1.html) | [Parte1.jl](Parte1.jl) |
| 2 | Campos vectoriales | [Abrir](Parte2.html) | [Parte2.jl](Parte2.jl) |
| 3 | Integrales de línea: campo vectorial | [Abrir](Parte3.html) | [Parte3.jl](Parte3.jl) |
| 4 | Teorema fundamental de las integrales de línea | [Abrir](Parte4.html) | [Parte4.jl](Parte4.jl) |
| 5 | Teorema de Green en el plano | [Abrir](Parte5.html) | [Parte5.jl](Parte5.jl) |
| 6 | Áreas de superficies | [Abrir](Parte6.html) | [Parte6.jl](Parte6.jl) |
| 7 | Integrales de superficie | [Abrir](Parte7.html) | [Parte7.jl](Parte7.jl) |
| 8 | Teorema de Stokes | [Abrir](Parte8.html) | [Parte8.jl](Parte8.jl) |
| 9 | Teorema de la divergencia | [Abrir](Parte9.html) | [Parte9.jl](Parte9.jl) |

## Visualizar las tareas

Abre [la página web](https://jostinbb.github.io/Calculo-Vectorial/) y pulsa **Abrir notebook** en la parte que quieras consultar. Los archivos HTML exportados por Pluto contienen los resultados guardados, las fórmulas y las gráficas. No requieren Julia para su lectura. Las nueve páginas tienen una barra flotante para regresar al índice y cambiar entre claro y oscuro, que se recuerda entre páginas.

Los archivos `.jl` contienen los notebooks originales y permiten volver a ejecutarlos desde Pluto.jl con las dependencias correspondientes. La web publicada es una **versión estática**, no un servidor de Julia.

## Publicar en GitHub Pages

1. Crea un repositorio público en GitHub, por ejemplo `Calculo-Vectorial-Taller-1`.
2. Sube **el contenido de esta carpeta**, no el ZIP completo ni una carpeta anidada. Debe quedar `index.html` visible en la raíz del repositorio.
3. En el repositorio, abre **Settings → Pages**.
4. En **Build and deployment → Source**, selecciona **Deploy from a branch**.
5. Selecciona **main** y **/(root)** y guarda.
6. Abre la web en [https://jostinbb.github.io/Calculo-Vectorial/](https://jostinbb.github.io/Calculo-Vectorial/).

La presencia de `.nojekyll` permite publicar los archivos estáticos directamente. No necesitas un servidor, dominio de pago ni una compilación adicional.

Más indicaciones en [COMO_PUBLICAR.md](COMO_PUBLICAR.md).

## Estructura

```text
Calculo-Vectorial/
├── index.html            # Portada y navegación (estilo Pluto)
├── favicon.svg           # Icono del sitio
├── pluto-site-theme.css  # Colores originales de Pluto, seleccionables
├── pluto-site-theme.js   # Alternancia de tema en las nueve partes
├── .nojekyll             # Sitio estático para GitHub Pages
├── README.md             # Información y enlaces
├── COMO_PUBLICAR.md      # Instrucciones de publicación
├── Parte1.html ... Parte9.html  # Exportaciones originales de Pluto
└── Parte1.jl   ... Parte9.jl    # Código Julia original
```

> **Nota:** el README original del archivo comprimido tenía temas mal asignados. Esta tabla se corrigió consultando los encabezados de cada `.jl`, sin modificar el contenido académico de los notebooks. Los HTML exportados recibieron únicamente el enlace a los estilos y controles de navegación.
