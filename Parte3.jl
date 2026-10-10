### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ b22e2d97-cffe-40e8-a857-a4430b51bdfd
using SymPy

# ╔═╡ 10772dcb-3400-4b3d-a03b-966d63ea5a48
using LinearAlgebra

# ╔═╡ 9be90baf-368c-4cd8-8f41-643154a959b6
@syms t::real x::real y::real z::real 

# ╔═╡ 14fddbb6-ed30-49b5-937c-6da4b237f5ee
md"# 3 Integrales de Línea: Campo Vectorial"

# ╔═╡ 6f06aaa1-d431-4ff6-b6a2-9a34f05393c0
md"## 3.1 Encontrar la integral de línea de $\vec{F}$ desde $(0, 0, 0)$ hasta $(1, 1, 1)$ para cada una de las siguientes trayectorias:
a) $C_1: \vec{r}(t) = t\hat{i} + t\hat{j} + t\hat{k}, \quad 0 \le t \le 1$  
b) $C_2: \vec{r}(t) = t\hat{i} + t^2\hat{j} + t^4\hat{k}, \quad 0 \le t \le 1$  
c) Trayectoria $C_1 \cup C_2$ que consiste de las rectas $(0, 0, 0) \to (1, 1, 0)$ y $(1, 1, 0) \to (1, 1, 1)$"

# ╔═╡ bbe882f6-d800-4aa5-963e-00e695b99964
md"### 3.1.1  
$\vec{F} = 3y\hat{i} + 2x\hat{j} + 4z\hat{k}$"

# ╔═╡ c20f28b2-9c39-4375-ac66-b112339de2af
md"a"

# ╔═╡ 06934f51-ef38-4d8e-9f3f-b434d3e26e11
let
	r(t)=[t,t,t]
	F=[3t,2t,4t]
	dr=diff.(r(t))
	integrate(dot(F,dr),(t,0,1))
end

# ╔═╡ e4ebed0b-ce69-4551-b836-08e1df60ab12
md"b"

# ╔═╡ e52c38b4-5083-4764-ad7f-d9034f6fa576
let
	r(t)=[t,t^2,t^4]
	F=[3t^2,2t,4t^4]
	dr=diff.(r(t))
	integrate(dot(F,dr),(t,0,1))
end

# ╔═╡ d7a0917e-50b7-4c79-8410-d2590c379c8e
md"c"

# ╔═╡ 8e9b06bd-67ec-4d63-80d6-4e90ae15ec2b
let
	r(t)=[t,t,0]
		F=[3t,2t,0]
	dr=diff.(r(t))

	r2(t)=[1,1,t+1]
		F2=[3,2,4(t+1)]
	dr2=diff.(r2(t))

	
	integrate(dot(F,dr),(t,0,1))+	integrate(dot(F2,dr2),(t,0,1))
	
end

# ╔═╡ 2efdc7d0-a331-4e92-91d7-09527fd20e3f
md"### 3.1.2  
$\vec{F} = (3x^2 - 3x)\hat{i} + 3z\hat{j} + \hat{k}$"

# ╔═╡ 973bed33-d6f4-4f9b-aba0-33e9d277c8f2
md"a"

# ╔═╡ a828aa48-a554-43fe-8a61-97ac0ad924c3
let
	r(t)=[t,t,t]
	F=[3t^2 - 3t ,3t,1]
	dr=diff.(r(t))
	integrate(dot(F,dr),(t,0,1))
end

# ╔═╡ f5ca0a39-dced-428e-8371-a637a02a9eb2
md"b"

# ╔═╡ 88f236ee-2872-4a8e-96bc-1a5bc7b8052e
let
	r(t)=[t,t^2,t^4]
		F=[3t^2 - 3t ,3t^4,1]
	dr=diff.(r(t))
	integrate(dot(F,dr),(t,0,1))
end

# ╔═╡ 609c2503-72dc-4cb2-866c-7f17602722e4
md"c"

# ╔═╡ 120b8efc-f90b-4fc4-9b16-909ee91ffff8
let
	r(t)=[t,t,0]
		F=[3t^2 - 3t ,0,1]
	dr=diff.(r(t))

	r2(t)=[1,1,t+1]
	F2=[3 - 3 ,t+1,1]
	dr2=diff.(r2(t))

	
	integrate(dot(F,dr),(t,0,1))+	integrate(dot(F2,dr2),(t,0,1))
	
end

# ╔═╡ 47793774-1995-46c6-8acc-2ec336a18808
md"## 3.2 Evaluar las siguientes integrales considerando la trayectoria $C: \vec{r}(t) = t\hat{i} - \hat{j} + t^2\hat{k}, \quad 0 \le t \le 1$:"

# ╔═╡ ed2c933a-cdac-4378-ad5e-3f6acb4e7d89
md"a) $\int_{C} (x + y - z) \, dx$  
b) $\int_{C} (x + y - z) \, dz$"

# ╔═╡ fd404c5d-21b5-45ae-b71a-ea52576dae3a
let
	integrate(t-1-t^2,(t,0,1))
end

# ╔═╡ b2b3a6ae-c7b6-4d71-a1df-154d6d354fa5
let
	integrate((t-1-t^2)*2t,(t,0,1))
end

# ╔═╡ db663421-334e-48e4-8841-510ccdc2bae7
md"## 3.3 Trabajo realizado por $\vec{F}$
Encontrar el trabajo realizado por $\vec{F}$ sobre la curva en la dirección incremental de $t$."

# ╔═╡ d9547775-ebc4-4e87-b785-e8f4ff55f66b
md"### 3.3.1  
$\vec{F} = xy\hat{i} + y\hat{j} - yz\hat{k}, \quad \vec{r}(t) = t\hat{i} + t^2\hat{j} + t\hat{k}, \quad 0 \le t \le 1$"

# ╔═╡ f877dfa2-d9af-432e-8fbe-aa54dde043b3
let
	r(t)=[t,t^2,t]
	F=[3*t*t^2,t^2,-t^2 *t]
	dr=diff.(r(t))
	integrate(dot(F,dr),(t,0,1))
end

# ╔═╡ 6a40af40-9d47-4349-9dfa-6c6c99792b0c
md"### 3.3.2  
$\vec{F} = z\hat{i} + x\hat{j} + y\hat{k}, \quad \vec{r}(t) = \sin(t)\hat{i} + \cos(t)\hat{j} + t\hat{k}, \quad 0 \le t \le 2\pi$"

# ╔═╡ eb3898c0-e9e3-45e6-a7ab-2a77040147bb
let
	r(t)=[sin(t),cos(t),t]
	F=[t,sin(t),cos(t)]
	dr=diff.(r(t))
	integrate(dot(F,dr),(t,0,2*PI))
end

# ╔═╡ a3f30924-0862-45fc-bf9c-6ea0062f43c0
md"## 3.4 Evaluar la integral considerando la curva $y = x^2$ desde $(-1, 1)$ hasta $(2, 4)$:
$\int_{C} xy \, dx + (x + y) \, dy$"

# ╔═╡ cf626227-0c39-4018-8f5c-d14cb7b32b8d
let
	integrate((t^3+2t*(t+t^2)),(t,-1,2))
end

# ╔═╡ 07811664-8025-43a1-a6d9-39426b04e8f3
md"## 3.5 Evaluar $\int_{C} \vec{F} \cdot \vec{T} \, ds$ para el campo $\vec{F} = x^2\hat{i} - y\hat{j}$ sobre la trayectoria $x = y^2$ desde $(4, 2)$ hasta $(1, -1)$."

# ╔═╡ 00f4e47a-1bb1-43c2-ab79-cfdd6d866dd9
integrate(dot([t^4,-t],[2t,1]),(t,2,-1))

# ╔═╡ 20a02ebf-8d35-4aa3-adf1-46ad56d498d0
md"## 3.6 Circulación y Flux
Encontrar la circulación y el flujo de los campos $\vec{F}_1 = x\hat{i} + y\hat{j}$ y $\vec{F}_2 = -y\hat{i} + x\hat{j}$ alrededor y a través de las siguientes curvas:
a) Círculo $\vec{r}(t) = \cos(t)\hat{i} + \sin(t)\hat{j}, \quad 0 \le t \le 2\pi$  
b) Elipse $\vec{r}(t) = \cos(t)\hat{i} + 4\sin(t)\hat{j}, \quad 0 \le t \le 2\pi$"

# ╔═╡ 12c678a7-ddec-483b-8fe7-e94237cef4c9
md"a"

# ╔═╡ af85ce43-2aed-4a7e-bace-c88c81bac38c
md"F1"

# ╔═╡ 7a1d5b02-f259-4d52-adc4-649e0dab9695
md"circualcion"

# ╔═╡ 201399b2-5d89-40df-80c4-cc8bfce89fcb
integrate(dot([cos(t),sin(t)],[-sin(t),cos(t)]),(t,0,2pi))

# ╔═╡ 5ca74270-b1e4-4624-ae1e-8a902d4a4264
md"flux"

# ╔═╡ cf0c2578-858f-4288-aa44-acdccc6c18a8
integrate(dot([cos(t),sin(t)],[cos(t),sin(t)]),(t,0,2PI)) 

# ╔═╡ 67cabf5e-6f2c-43e0-8c2c-31678359c527
md"F2"

# ╔═╡ c6981786-6728-44cd-9446-c83d7c6b4ace
md"ciruclacion"

# ╔═╡ cf263ebe-79a1-407f-850b-9bc5a1b32f9f
integrate(dot([-sin(t),cos(t)],[-sin(t),cos(t)]),(t,0,2PI))

# ╔═╡ 88dd26b3-118f-420e-bba4-3ccf5c3631e5
md"flux"

# ╔═╡ c61e3c3b-ab41-4b6f-9661-b56f41db3c37
integrate(dot([-sin(t),cos(t)],[cos(t),sin(t)]),(t,0,2PI))

# ╔═╡ f6d7d1d1-e49c-4aab-97e8-0006aa9a371f
md"b"

# ╔═╡ 07483ef3-82e6-43cf-8923-c7266fc40d3a
md"F1"

# ╔═╡ 826a6762-57a2-42f2-9a9f-7244ad2247ad
md"circulacion"

# ╔═╡ 143a014e-348c-418d-8f71-e513c02720b7
integrate(dot([cos(t),4*sin(t)],[-sin(t),4*cos(t)]),(t,0,2PI))

# ╔═╡ db0a615f-827d-4036-9e09-6eff7f6088c3
md"flux"

# ╔═╡ 30ec9d32-f705-4149-9c62-084308118274
integrate(dot([cos(t),4*sin(t)],[4*cos(t),sin(t)]),(t,0,2PI))

# ╔═╡ c0144aee-eeda-4310-b3cf-8c17bbc0a359
md"F2"

# ╔═╡ 885f9a2f-5879-4ba0-96bf-dc4e0565bbdd
md"circulacion"

# ╔═╡ 4b16393e-fe7c-4df1-ab28-83a2ea0b22e0
integrate(dot([-4*sin(t),cos(t)],[-sin(t),4*cos(t)]),(t,0,2PI))

# ╔═╡ c154c055-45a3-4140-8f73-d1497db7017f
md"flux"

# ╔═╡ 2634e977-5bbe-4f7a-834d-6b1cf19e3bb2
integrate(dot([-4*sin(t),cos(t)],[4*cos(t),sin(t)]),(t,0,2PI))

# ╔═╡ 4f759f30-6b94-44a9-85fb-c6cc395228e6
md"## 3.7 Flujo del campo de velocidad $\vec{F} = y^2\hat{i} + 2xy\hat{j}$ sobre las trayectorias desde $(0, 0)$ hasta $(2, 4)$:
* Trayectoria 1: $y = x^2$  
* Trayectoria 2: $y = 2x$"

# ╔═╡ acf9be4f-5f71-488c-82b5-6341fbee20e5
integrate(dot([t^4,2*t^3],[2t,-1]),(t,0,2))

# ╔═╡ 0ac046b7-eed3-4f75-a3ab-106af45ea42e
integrate(dot([(2*t)^2,(2*t*2t)],[2,-1]),(t,0,2))

# ╔═╡ 559a70ad-4a32-4bc2-b9f0-511de9fffbb2
md"## 3.8 Flujo del campo de velocidad $\vec{F} = \frac{x}{y+1}\hat{i} + \frac{y}{x+1}\hat{j}$ sobre la curva $\vec{r}(t) = t^2\hat{i} + t\hat{j}, \quad 0 \le t \le 1$."

# ╔═╡ 79e1e90f-0d9f-4569-beca-7f082f1305b4
integrate(dot([t^2 /(t+1),(t/(t^2 +1))],[1,-2t]),(t,0,1))

# ╔═╡ 608478b5-60ae-4e7c-93ba-67e3e53dcfa9
md"## 3.9 Campo de vectores unitarios
Encontrar un campo $\vec{F} = M(x, y)\hat{i} + N(x, y)\hat{j}$ en el plano $xy$ tal que para cada punto $(x, y) \neq (0, 0)$, $\vec{F}$ sea un vector unitario apuntando hacia el origen."

# ╔═╡ e6661cec-d0f1-4619-89e9-aedfd6196281
let
	F=[-x,-y]
end

# ╔═╡ 38f2ac34-7d2d-407d-81d0-ddf47e0af492
md"## 3.10 Flujo de fluido en el espacio
$\vec{F}$ es el campo de velocidad. Encontrar el flujo sobre la curva definida en la dirección incremental de $t$."

# ╔═╡ 6037d3db-7171-42a2-8675-a3fe7b6cb0f7
md"### 3.10.1  
$\vec{F} = -4xy\hat{i} + 8y\hat{j} + 2\hat{k}, \quad \vec{r}(t) = t\hat{i} + t^2\hat{j} + \hat{k}, \quad 0 \le t \le 2$"

# ╔═╡ 47ce4772-a275-4817-8f83-393a88e59a0e
let
	F(x,y)=[-4x*y,8y,2]
	r=[t,t^2,1]
	dr=diff.(r)
	integrate(dot(F(t,t^2),dr),(t,0,2))
end

# ╔═╡ e33f9fa8-f795-48cf-bdfb-09e37b72e748
md"### 3.10.2  
$\vec{F} = (x - z)\hat{i} + x\hat{k}, \quad \vec{r}(t) = \cos(t)\hat{i} + \sin(t)\hat{k}, \quad 0 \le t \le \pi$"

# ╔═╡ d4e88779-91d6-486c-aea0-7af6679a4d5a
let
	F(x,z)=[x-z,0,x]
	r=[cos(t),0,sin(t)]
	dr=diff.(r)
	integrate(dot(F(cos(t),sin(t)),dr),(t,0,PI))
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[compat]
SymPy = "~2.3.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.6"
manifest_format = "2.0"
project_hash = "a741edd54e859e0f84e024e22ac3afb6bbfcda6d"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.CommonEq]]
git-tree-sha1 = "6b0f0354b8eb954cdba708fb262ef00ee7274468"
uuid = "3709ef60-1bee-4518-9f2f-acd86f176c50"
version = "0.2.1"

[[deps.CommonSolve]]
git-tree-sha1 = "78ea4ddbcf9c241827e7035c3a03e2e456711470"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.6"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "8f06b0cfa4c514c7b9546756dbae91fcfbc92dc9"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.10.3"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "b3ad4a0255688dcb895a52fafbaae3023b588a90"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.4.0"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6893345fd6658c8e475d40155789f4860ac3b21"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.4+0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.11.4"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "07a921781cab75691315adc645096ed5e370cb77"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "9816a3826b0ebf49ab4926e2b18842ad8b5c8f04"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.96.4"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "5acc6a41b3082920f79ca3c759acbcecf18a8d78"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.7.1"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "fa95b3b097bcef5845c142ea2e085f1b2591e92c"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.7.1"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsStaticArraysCoreExt = ["StaticArraysCore"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SymPy]]
deps = ["CommonEq", "CommonSolve", "LinearAlgebra", "PyCall", "SpecialFunctions", "SymPyCore"]
git-tree-sha1 = "d3c2de8adc6e36352d2a2dae3ae87099964fcbc0"
uuid = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
version = "2.3.3"

[[deps.SymPyCore]]
deps = ["CommonEq", "CommonSolve", "Latexify", "LinearAlgebra", "Markdown", "RecipesBase", "SpecialFunctions", "TermInterface"]
git-tree-sha1 = "efca83098a5badccda86498fd692a88d26c12ade"
uuid = "458b697b-88f0-4a86-b56b-78b75cfb3531"
version = "0.3.3"

    [deps.SymPyCore.extensions]
    SymPyCoreSymbolicUtilsExt = "SymbolicUtils"

    [deps.SymPyCore.weakdeps]
    SymbolicUtils = "d1185830-fcd6-423d-90d6-eec64667417b"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TermInterface]]
git-tree-sha1 = "d673e0aca9e46a2f63720201f55cc7b3e7169b16"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "2.0.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"
"""

# ╔═╡ Cell order:
# ╠═b22e2d97-cffe-40e8-a857-a4430b51bdfd
# ╠═10772dcb-3400-4b3d-a03b-966d63ea5a48
# ╠═9be90baf-368c-4cd8-8f41-643154a959b6
# ╟─14fddbb6-ed30-49b5-937c-6da4b237f5ee
# ╟─6f06aaa1-d431-4ff6-b6a2-9a34f05393c0
# ╟─bbe882f6-d800-4aa5-963e-00e695b99964
# ╟─c20f28b2-9c39-4375-ac66-b112339de2af
# ╠═06934f51-ef38-4d8e-9f3f-b434d3e26e11
# ╟─e4ebed0b-ce69-4551-b836-08e1df60ab12
# ╠═e52c38b4-5083-4764-ad7f-d9034f6fa576
# ╟─d7a0917e-50b7-4c79-8410-d2590c379c8e
# ╠═8e9b06bd-67ec-4d63-80d6-4e90ae15ec2b
# ╟─2efdc7d0-a331-4e92-91d7-09527fd20e3f
# ╟─973bed33-d6f4-4f9b-aba0-33e9d277c8f2
# ╠═a828aa48-a554-43fe-8a61-97ac0ad924c3
# ╟─f5ca0a39-dced-428e-8371-a637a02a9eb2
# ╠═88f236ee-2872-4a8e-96bc-1a5bc7b8052e
# ╟─609c2503-72dc-4cb2-866c-7f17602722e4
# ╠═120b8efc-f90b-4fc4-9b16-909ee91ffff8
# ╟─47793774-1995-46c6-8acc-2ec336a18808
# ╟─ed2c933a-cdac-4378-ad5e-3f6acb4e7d89
# ╠═fd404c5d-21b5-45ae-b71a-ea52576dae3a
# ╠═b2b3a6ae-c7b6-4d71-a1df-154d6d354fa5
# ╟─db663421-334e-48e4-8841-510ccdc2bae7
# ╟─d9547775-ebc4-4e87-b785-e8f4ff55f66b
# ╠═f877dfa2-d9af-432e-8fbe-aa54dde043b3
# ╟─6a40af40-9d47-4349-9dfa-6c6c99792b0c
# ╠═eb3898c0-e9e3-45e6-a7ab-2a77040147bb
# ╟─a3f30924-0862-45fc-bf9c-6ea0062f43c0
# ╠═cf626227-0c39-4018-8f5c-d14cb7b32b8d
# ╟─07811664-8025-43a1-a6d9-39426b04e8f3
# ╠═00f4e47a-1bb1-43c2-ab79-cfdd6d866dd9
# ╟─20a02ebf-8d35-4aa3-adf1-46ad56d498d0
# ╟─12c678a7-ddec-483b-8fe7-e94237cef4c9
# ╟─af85ce43-2aed-4a7e-bace-c88c81bac38c
# ╟─7a1d5b02-f259-4d52-adc4-649e0dab9695
# ╠═201399b2-5d89-40df-80c4-cc8bfce89fcb
# ╟─5ca74270-b1e4-4624-ae1e-8a902d4a4264
# ╠═cf0c2578-858f-4288-aa44-acdccc6c18a8
# ╟─67cabf5e-6f2c-43e0-8c2c-31678359c527
# ╟─c6981786-6728-44cd-9446-c83d7c6b4ace
# ╠═cf263ebe-79a1-407f-850b-9bc5a1b32f9f
# ╟─88dd26b3-118f-420e-bba4-3ccf5c3631e5
# ╠═c61e3c3b-ab41-4b6f-9661-b56f41db3c37
# ╟─f6d7d1d1-e49c-4aab-97e8-0006aa9a371f
# ╟─07483ef3-82e6-43cf-8923-c7266fc40d3a
# ╟─826a6762-57a2-42f2-9a9f-7244ad2247ad
# ╠═143a014e-348c-418d-8f71-e513c02720b7
# ╟─db0a615f-827d-4036-9e09-6eff7f6088c3
# ╠═30ec9d32-f705-4149-9c62-084308118274
# ╟─c0144aee-eeda-4310-b3cf-8c17bbc0a359
# ╟─885f9a2f-5879-4ba0-96bf-dc4e0565bbdd
# ╠═4b16393e-fe7c-4df1-ab28-83a2ea0b22e0
# ╟─c154c055-45a3-4140-8f73-d1497db7017f
# ╠═2634e977-5bbe-4f7a-834d-6b1cf19e3bb2
# ╟─4f759f30-6b94-44a9-85fb-c6cc395228e6
# ╠═acf9be4f-5f71-488c-82b5-6341fbee20e5
# ╠═0ac046b7-eed3-4f75-a3ab-106af45ea42e
# ╟─559a70ad-4a32-4bc2-b9f0-511de9fffbb2
# ╠═79e1e90f-0d9f-4569-beca-7f082f1305b4
# ╟─608478b5-60ae-4e7c-93ba-67e3e53dcfa9
# ╠═e6661cec-d0f1-4619-89e9-aedfd6196281
# ╟─38f2ac34-7d2d-407d-81d0-ddf47e0af492
# ╟─6037d3db-7171-42a2-8675-a3fe7b6cb0f7
# ╠═47ce4772-a275-4817-8f83-393a88e59a0e
# ╟─e33f9fa8-f795-48cf-bdfb-09e37b72e748
# ╠═d4e88779-91d6-486c-aea0-7af6679a4d5a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
