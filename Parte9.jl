### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ e8bb22fa-63fa-4881-bab9-544620ff2bbb
using SymPy

# ╔═╡ af94cb9d-1260-456e-85fe-15113e6e9074
using LinearAlgebra

# ╔═╡ ffebfae8-aec5-4e70-94f3-95bbd8f705bc
md"# 9 Teorema de la Divergencia"

# ╔═╡ 0fd3c13b-9e5f-4a99-9398-18ac2a574ff5
@syms x::Real y::Real z::Real t::Real θ::Real r::Real ϕ::Real u::Real v::Real c::Real a::Real b::Real p::Real

# ╔═╡ 46d42891-6668-4601-8c5c-3c621a3204a2
md"## 9.1 Divergencia de Campos Vectoriales
Encontrar la divergencia $\nabla \cdot \vec{F}$ de cada campo vectorial:"

# ╔═╡ 75626c3e-095a-42c8-a81e-46eaeeeb77ea
md"### 9.1.1  
$\vec{F} = (x - y + z)\hat{i} + (2x + y - z)\hat{j} + (3x + 2y - 2z)\hat{k}$"

# ╔═╡ 9a8e863b-be63-44f4-8545-3ebda591ab74
let
	F=[x-y+z,2x+y-z,3x+2y-2z]
	divF=diff(F[1],x)+diff(F[2],y)+diff(F[3],z)
end	

# ╔═╡ 16579199-50fe-4579-a5b1-de1cf950a3ba
md"### 9.1.2  
$\vec{F} = ye^{xyz}\hat{i} + ze^{xyz}\hat{j} + xe^{xyz}\hat{k}$"

# ╔═╡ 35bcfcaf-628b-4570-bba3-edadb259402e
let
	F = [y*exp(x*y*z), z*exp(x*y*z), x*exp(x*y*z)]
	divF = diff(F[1],x) + diff(F[2],y) + diff(F[3],z)
end

# ╔═╡ e8d249ff-4eee-4ee8-bcb3-0f265d5b6447
md"### 9.1.3  
$\vec{F} = \frac{-y\hat{i} + x\hat{j}}{(x^2 + y^2)^{1/2}}$"

# ╔═╡ 76765b18-829e-41d8-a726-dff7e1a9a01b
let

	F = [-y/(x^2 + y^2)^(1/2), x/(x^2 + y^2)^(1/2)]
	divF = diff(F[1],x) + diff(F[2],y) 
end

# ╔═╡ f366ed36-dcd0-440a-8532-89a1e6b88a1d
md"## 9.2 Cálculo de Flux mediante el Teorema de la Divergencia
Calcular el Flux hacia afuera de $\vec{F}$ a través de la región de frontera $D$ utilizando el Teorema de la Divergencia:
$$\iint_{S} \vec{F} \cdot \vec{n} \, d\sigma = \iiint_{D} (\nabla \cdot \vec{F}) \, dV$$"

# ╔═╡ 8b952ef7-a166-4b33-9cfe-0786b74046d5
md"### 9.2.1  
$\vec{F} = (y - x)\hat{i} + (z - y)\hat{j} + (y - x)\hat{k}$. $D$: Cubo encerrado por los planos $x = \pm 1, y = \pm 1, z = \pm 1$."

# ╔═╡ 515e648a-18e1-449f-adc3-668fa8d8632b
let
	F=[y-x,z-y,y-x]
	divF=diff(F[1],x)+diff(F[2],y)+diff(F[3],z)
	integrate(divF,(x,-1,1),(y,-1,1),(z,-1,1))
end	

# ╔═╡ 9f632d51-6ce0-4e98-93bc-a1bc80ec0863
md"### 9.2.2  
$\vec{F} = y\hat{i} + xy\hat{j} - z\hat{k}$. $D$: Región dentro del cilindro sólido $x^2 + y^2 \le 4$, entre los planos $z = 0$ y el paraboloide $z = x^2 + y^2$."

# ╔═╡ 60aba21a-f441-4628-a3b7-35a6ea99952a
let
	F=[y,x*y,-z]  
	divF=diff(F[1],x)+diff(F[2],y)+diff(F[3],z)
	integrate((r*cos(t)-1)*r,(r,0,2),(z,0,2),(t,0,2PI))
end	

# ╔═╡ 1e32ab31-a7f9-4524-b7e5-a25de66c4ffb
md"### 9.2.3  
$\vec{F} = \sqrt{x^2 + y^2 + z^2} (x\hat{i} + y\hat{j} + z\hat{k})$. $D$: Región esférica $1 \le x^2 + y^2 + z^2 \le 2$."

# ╔═╡ 4d720a1e-add4-4f4b-aff5-e56a8562c8f5
let
	F= sqrt(x^2+y^2+z^2) .*[x,y,z] 
		divF=diff(F[1],x)+diff(F[2],y)+diff(F[3],z)
	
end	

# ╔═╡ 9b9a09bb-72c9-4a00-873a-1610dea1a0e4
let
	
	integrate(((sin(p)*cos(t))^2+(sin(p)*sin(t))^2+(cos(p))^2+3r)*r^2*sin(p),(r,1,sqrt(2)),(p,0,PI),(t,0,2PI))
end

# ╔═╡ da82da44-0640-4f4c-a22d-48446f9b970d
md"### 9.2.4  
$\vec{F} = (5x^3 + 12xy^2)\hat{i} + (y^3 + e^y \sin z)\hat{j} + (5z^3 + e^y \cos z)\hat{k}$. $D$: Región sólida entre las esferas $x^2 + y^2 + z^2 = 1$ y $x^2 + y^2 + z^2 = 2$."

# ╔═╡ ba0f664a-5894-4ede-9c09-c77514e134a8
let
	F=[5x^3+12*x*y^2,y^3+exp(y)*sin(z),5z^3+exp(y)*cos(z)]
			divF=diff(F[1],x)+diff(F[2],y)+diff(F[3],z)

	integrate(15r^4*sin(p),(r,1,sqrt(2)),(p,0,PI),(t,0,2PI))
end

# ╔═╡ 52aec9e0-9c5f-4cbb-a5ee-2c412e243eae
md"## 9.3 Existencia de Potencial Vectorial
Si $\vec{F} = (y \cos(2x))\hat{i} + (y^2 \sin(2x))\hat{j} + (x^2 y + z)\hat{k}$, ¿existe un campo vectorial $\vec{A}$ tal que $\vec{F} = \nabla \times \vec{A}$? Explica tu respuesta basándote en las propiedades de la divergencia."

# ╔═╡ 49cbaa80-6529-4afc-97b8-8a9888ce4727
let
	F=[y*cos(2x),y^2*sin(2x),(x^2*y+z)]
			divF=diff(F[1],x)+diff(F[2],y)+diff(F[3],z)
end

# ╔═╡ 741e0af8-9afb-47d3-8d08-986739f05d56
md"Como la divergencia nos dio $1$ (y no $0$), NO existe un campo vectorial $\vec{A}$ tal que $\vec{F} = \nabla \times \vec{A}$"

# ╔═╡ 60b26597-2ad0-4e3a-abfa-18d4dbd9c7f7
md"## 9.4 Flux Neto en Regiones Complejas
Calcular el flux neto hacia el exterior del campo vectorial $\vec{F}$ sobre la superficie $S$ que rodea la región $D$ delimitada por los planos $y = 0, z = 0, z = 2-y$ y el cilindro parabólico $z = 1 - x^2$.
$$\vec{F} = xy\hat{i} + (\sin(xz) + y^2)\hat{j} + (e^{xy^2} + x)\hat{k}$$"

# ╔═╡ d04796f7-5cac-4bf1-b0e9-f3bb491f01e3
let
	F = [x*y, sin(x*z) + y^2, exp(x*y^2) + x]
	divF = diff(F[1], x) + diff(F[2], y) + diff(F[3], z)
	
	
	integrate(divF, (y, 0, 2-z), (z, 0, 1-x^2), (x, -1, 1))
end

# ╔═╡ e22cd6ca-85fd-43df-8d46-567b65f9c32d
md"## 9.5 Identidades Vectoriales
Si $\vec{F}$ es un campo vectorial diferenciable y $g(x, y, z)$ es una función escalar diferenciable, verificar las siguientes identidades:
a) $\nabla \cdot (g\vec{F}) = g(\nabla \cdot \vec{F}) + \nabla g \cdot \vec{F}$  
b) $\nabla \times (g\vec{F}) = g(\nabla \times \vec{F}) + \nabla g \times \vec{F}$"

# ╔═╡ fd90d51e-165e-478d-a3c1-f61d0cb78a41
g = SymFunction("g")

# ╔═╡ 74cbc093-ec83-4c66-abbe-56ce74168fc0
F1 = SymFunction("F1")

# ╔═╡ fe62e1d1-4d38-4fe1-b131-179d546a23d6
F2 = SymFunction("F2")

# ╔═╡ 309ce657-21a4-42fe-8d7c-9ae20666d5ed
F3 = SymFunction("F3")

# ╔═╡ a3e749f8-d4e7-42ab-9f00-7f0e4e5e2d09
g_xyz = g(x, y, z)

# ╔═╡ afe17dad-c3dd-4fef-8a2e-26f9b5aa4b93
F = [F1(x, y, z), F2(x, y, z), F3(x, y, z)]

# ╔═╡ 7e3b6324-546f-4e20-8f00-5d3791bd2b2c
md"Identidad de la Divergencia"

# ╔═╡ 3187c30d-a356-43a2-af35-88dc3562307f
gF = g_xyz .* F

# ╔═╡ 23e6fd22-bd76-4980-9db9-08914de51735
LHS_a = diff(gF[1], x) + diff(gF[2], y) + diff(gF[3], z)

# ╔═╡ 5727e692-c75d-4732-a441-8264b3497d81
divF = diff(F[1], x) + diff(F[2], y) + diff(F[3], z)

# ╔═╡ f9638b87-c7ef-470e-b4aa-4d0bafa669f3
grad_g = [diff(g_xyz, x), diff(g_xyz, y), diff(g_xyz, z)]

# ╔═╡ 4159f86e-1060-40e7-90b6-66e03b8b1679
RHS_a = g_xyz * divF + (grad_g[1]*F[1] + grad_g[2]*F[2] + grad_g[3]*F[3])

# ╔═╡ e1563898-0700-42ff-a7d1-ec27929b25d5
md" Comprobamos la igualdad restando y simplificando"

# ╔═╡ e6da0786-2b0b-4ca2-a1f5-ab94a37f8179
identidad_a = simplify(LHS_a - RHS_a) == 0

# ╔═╡ 510fe587-55ed-4ed2-aa13-348cfd907834
md" Identidad del Rotacional "

# ╔═╡ 4ceeb140-0e04-4094-adac-1affc772f20c
curl_LHS = [
		diff(gF[3], y) - diff(gF[2], z),
		diff(gF[1], z) - diff(gF[3], x),
		diff(gF[2], x) - diff(gF[1], y)
	]

# ╔═╡ 2c8707e7-27db-4e03-bbe2-b5c80fc3f099
curl_F = [
		diff(F[3], y) - diff(F[2], z),
		diff(F[1], z) - diff(F[3], x),
		diff(F[2], x) - diff(F[1], y)
	]

# ╔═╡ ba5e2bd8-2d9f-482c-81f4-793fa777deaa
cross_grad_F = [
		grad_g[2]*F[3] - grad_g[3]*F[2],
		grad_g[3]*F[1] - grad_g[1]*F[3],
		grad_g[1]*F[2] - grad_g[2]*F[1]
	]

# ╔═╡ 50f9b80b-4e1a-43ba-803b-0cdf0b8bf29f
curl_RHS = g_xyz .* curl_F .+ cross_grad_F

# ╔═╡ 8fc10dee-aabc-4c67-831e-1e2b9a86410d
md" Comprobamos la igualdad vector a vector"

# ╔═╡ b629ef8f-ab1d-4690-b3c6-79ef8a892a29
identidad_b = simplify.(curl_LHS .- curl_RHS) == [0, 0, 0]

# ╔═╡ 3f9b9638-69d5-444b-88d2-453d4f1731db
md"Retornamos el resultado booleano de ambas demostraciones"

# ╔═╡ 1137e965-66e5-4330-b0d4-40c1375a4a57
(Identidad_A_Verificada = identidad_a, Identidad_B_Verificada = identidad_b)

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
# ╟─ffebfae8-aec5-4e70-94f3-95bbd8f705bc
# ╠═e8bb22fa-63fa-4881-bab9-544620ff2bbb
# ╠═af94cb9d-1260-456e-85fe-15113e6e9074
# ╠═0fd3c13b-9e5f-4a99-9398-18ac2a574ff5
# ╟─46d42891-6668-4601-8c5c-3c621a3204a2
# ╟─75626c3e-095a-42c8-a81e-46eaeeeb77ea
# ╠═9a8e863b-be63-44f4-8545-3ebda591ab74
# ╟─16579199-50fe-4579-a5b1-de1cf950a3ba
# ╠═35bcfcaf-628b-4570-bba3-edadb259402e
# ╟─e8d249ff-4eee-4ee8-bcb3-0f265d5b6447
# ╠═76765b18-829e-41d8-a726-dff7e1a9a01b
# ╟─f366ed36-dcd0-440a-8532-89a1e6b88a1d
# ╟─8b952ef7-a166-4b33-9cfe-0786b74046d5
# ╠═515e648a-18e1-449f-adc3-668fa8d8632b
# ╟─9f632d51-6ce0-4e98-93bc-a1bc80ec0863
# ╠═60aba21a-f441-4628-a3b7-35a6ea99952a
# ╟─1e32ab31-a7f9-4524-b7e5-a25de66c4ffb
# ╠═4d720a1e-add4-4f4b-aff5-e56a8562c8f5
# ╠═9b9a09bb-72c9-4a00-873a-1610dea1a0e4
# ╟─da82da44-0640-4f4c-a22d-48446f9b970d
# ╠═ba0f664a-5894-4ede-9c09-c77514e134a8
# ╟─52aec9e0-9c5f-4cbb-a5ee-2c412e243eae
# ╠═49cbaa80-6529-4afc-97b8-8a9888ce4727
# ╟─741e0af8-9afb-47d3-8d08-986739f05d56
# ╟─60b26597-2ad0-4e3a-abfa-18d4dbd9c7f7
# ╠═d04796f7-5cac-4bf1-b0e9-f3bb491f01e3
# ╟─e22cd6ca-85fd-43df-8d46-567b65f9c32d
# ╠═fd90d51e-165e-478d-a3c1-f61d0cb78a41
# ╠═74cbc093-ec83-4c66-abbe-56ce74168fc0
# ╠═fe62e1d1-4d38-4fe1-b131-179d546a23d6
# ╠═309ce657-21a4-42fe-8d7c-9ae20666d5ed
# ╠═a3e749f8-d4e7-42ab-9f00-7f0e4e5e2d09
# ╠═afe17dad-c3dd-4fef-8a2e-26f9b5aa4b93
# ╠═7e3b6324-546f-4e20-8f00-5d3791bd2b2c
# ╠═3187c30d-a356-43a2-af35-88dc3562307f
# ╠═23e6fd22-bd76-4980-9db9-08914de51735
# ╠═5727e692-c75d-4732-a441-8264b3497d81
# ╠═f9638b87-c7ef-470e-b4aa-4d0bafa669f3
# ╠═4159f86e-1060-40e7-90b6-66e03b8b1679
# ╠═e1563898-0700-42ff-a7d1-ec27929b25d5
# ╠═e6da0786-2b0b-4ca2-a1f5-ab94a37f8179
# ╠═510fe587-55ed-4ed2-aa13-348cfd907834
# ╠═4ceeb140-0e04-4094-adac-1affc772f20c
# ╠═2c8707e7-27db-4e03-bbe2-b5c80fc3f099
# ╠═ba5e2bd8-2d9f-482c-81f4-793fa777deaa
# ╠═50f9b80b-4e1a-43ba-803b-0cdf0b8bf29f
# ╠═8fc10dee-aabc-4c67-831e-1e2b9a86410d
# ╠═b629ef8f-ab1d-4690-b3c6-79ef8a892a29
# ╠═3f9b9638-69d5-444b-88d2-453d4f1731db
# ╠═1137e965-66e5-4330-b0d4-40c1375a4a57
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
