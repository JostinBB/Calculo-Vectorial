### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ bcfe4bdc-6615-4d82-8b3f-ec187f4ca456
md"# 4 Teorema Fundamental de las Integrales de Línea"

# ╔═╡ 41710abc-0d51-443d-a3a8-27860348837f
md"## 4.1 Campos Vectoriales Conservativos
Determinar si $\vec{F}$ es un campo vectorial conservativo. Si es así, encontrar una función potencial $f$ tal que $\vec{F} = \nabla f$."

# ╔═╡ 45bea9e5-2c7c-4f88-83c3-e65b3c35d7d5
md"### 4.1.1  
$\vec{F}(x, y) = (2x - 3y)\hat{i} + (-3x + 4y - 8)\hat{j}$"

# ╔═╡ ff49ddb7-2147-4e4e-b63f-1a2b456e8f66
md"### 4.1.2  
$\vec{F}(x, y) = (ye^{x} + \sin(y))\hat{i} + (e^{x} + x\cos(y))\hat{j}$"

# ╔═╡ 0a6c04da-e9b2-4f35-81d2-91a58d296957
md"## 4.2 Independencia de la Trayectoria
Demostrar que la integral de línea es independiente de la trayectoria y evaluar la integral:
$$\int_{C} 2xe^{-y} \, dx + (2y - x^2 e^{-y}) \, dy$$
donde $C$ es cualquier trayectoria desde $(1, 0)$ hasta $(2, 1)$."

# ╔═╡ e82550e9-dc75-4958-9413-350f4b3bace4
md"## 4.3 Propiedades de Campos Conservativos en $\mathbb{R}^3$
Demostrar que si el campo vectorial $\vec{F} = P\hat{i} + Q\hat{j} + R\hat{k}$ es conservativo y $P, Q, R$ tienen derivadas parciales continuas de primer orden, entonces:
$$\frac{\partial P}{\partial y} = \frac{\partial Q}{\partial x}, \quad \frac{\partial P}{\partial z} = \frac{\partial R}{\partial x}, \quad \frac{\partial Q}{\partial z} = \frac{\partial R}{\partial y}$$"

# ╔═╡ 9d9e2abb-b420-4319-9560-6e55c6bf4174
md"## 4.4 Topología de Conjuntos
Determinar si el conjunto dado es: **a) abierto**, **b) conexo** y **c) simplemente conexo**."

# ╔═╡ 9baece73-9d74-4e59-b17b-f600f9c70a56
md"### 4.4.1  
$\{(x, y) \mid 0 < y < 3\}$"

# ╔═╡ 66669784-90a6-433b-88d1-932d225d7e2d
md"### 4.4.2  
$\{(x, y) \mid 1 \le x^2 + y^2 \le 4, \ y \ge 0\}$"

# ╔═╡ bae97f2d-5fd3-4320-8f2d-1531d4137b04
md"### 4.4.3  
$\{(x, y) \mid (x, y) \neq (2, 3)\}$"

# ╔═╡ 011cb45d-5ce4-4962-889b-2a050f764e2d
md"## 4.5 Análisis de un Campo No Conservativo Especial
Sea $\vec{F}(x, y) = P(x, y)\hat{i} + Q(x, y)\hat{j} = \frac{-y\hat{i} + x\hat{j}}{x^2 + y^2}$

a) Demostrar que $\frac{\partial P}{\partial y} = \frac{\partial Q}{\partial x}$.  
b) Demostrar que $\int_{C} \vec{F} \cdot d\vec{r}$ **no** es independiente de la trayectoria.  
*Sugerencia: calcular las integrales sobre $C_1$ (mitad superior) y $C_2$ (mitad inferior) de la circunferencia $x^2 + y^2 = 1$ desde $(1, 0)$ hasta $(-1, 0)$.*"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.5"
manifest_format = "2.0"
project_hash = "71853c6197a6a7f222db0f1978c7cb232b87c5ee"

[deps]
"""

# ╔═╡ Cell order:
# ╠═bcfe4bdc-6615-4d82-8b3f-ec187f4ca456
# ╠═41710abc-0d51-443d-a3a8-27860348837f
# ╠═45bea9e5-2c7c-4f88-83c3-e65b3c35d7d5
# ╠═ff49ddb7-2147-4e4e-b63f-1a2b456e8f66
# ╠═0a6c04da-e9b2-4f35-81d2-91a58d296957
# ╠═e82550e9-dc75-4958-9413-350f4b3bace4
# ╠═9d9e2abb-b420-4319-9560-6e55c6bf4174
# ╠═9baece73-9d74-4e59-b17b-f600f9c70a56
# ╠═66669784-90a6-433b-88d1-932d225d7e2d
# ╠═bae97f2d-5fd3-4320-8f2d-1531d4137b04
# ╠═011cb45d-5ce4-4962-889b-2a050f764e2d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
