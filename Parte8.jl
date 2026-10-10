### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ c6f01dbc-9a83-492e-9389-16776fb76c3d
using SymPy

# ╔═╡ 7c210c25-07b9-400b-b75e-bd4c9c8d5e68
using LinearAlgebra

# ╔═╡ 2be2dcd5-987b-4677-86a6-782f9ce6ff02
@syms x::Real y::Real z::Real t::Real θ::Real r::Real ϕ::Real u::Real v::Real c::Real a::Real b::Real p::Real

# ╔═╡ d34d4f7a-1832-4b18-a6d0-01c41063146a
md"# 8 Teorema de Stokes"

# ╔═╡ eedce593-488f-4707-870c-ebc33ea96641
md"## 8.1 Rotacional de Campos Vectoriales
Encontrar el curl (rotacional) $\nabla \times \vec{F}$ para cada campo vectorial:"

# ╔═╡ 001176ce-cff7-49ff-b1d7-47a40885faf8
md"### 8.1.1  
$\vec{F} = (x + y - z)\hat{i} + (2x - y + 3z)\hat{j} + (3x + 2y + z)\hat{k}$"

# ╔═╡ ea8109f3-7931-42be-bf42-cc994fdcf62c
let
	F= [x+y-z,2*x-y+3z,3x+2y+z]
	rotF=[
 diff(F[3], y) - diff(F[2], z),
 diff(F[1], z) - diff(F[3], x),
 diff(F[2], x) - diff(F[1], y)
]
	
end

# ╔═╡ f1ecd628-a2a8-4daa-bd84-e6435f320935
md"### 8.1.2  
$\vec{F} = (xy + z)\hat{i} + (yz + x)\hat{j} + (xz + y)\hat{k}$"

# ╔═╡ d04316c5-964a-478e-b9b4-6e8f64cf8115
let
	F= [x*y+z,y*z+x,x*z+y]
	rotF=[
 diff(F[3], y) - diff(F[2], z),
 diff(F[1], z) - diff(F[3], x),
 diff(F[2], x) - diff(F[1], y)
]
end

# ╔═╡ 42c94328-3d8b-486a-9fed-8378171904dc
md"### 8.1.3  
$\vec{F} = x^2yz\hat{i} + xy^2z\hat{j} + xyz^2\hat{k}$"

# ╔═╡ 4d8c93ba-a26c-488f-a221-57e831a02daf
let
	F= [x^2*y*z,x*y^2*z,x*y*z^2]
rotF=[
 diff(F[3], y) - diff(F[2], z),
 diff(F[1], z) - diff(F[3], x),
 diff(F[2], x) - diff(F[1], y)
]
	
end

# ╔═╡ a23a72c8-5862-45e4-8512-19545b21bd48
md"## 8.2 Circulación mediante el Teorema de Stokes
Emplear la integral de superficie del Teorema de Stokes para calcular la circulación de $\vec{F}$ sobre la curva $C$ en la dirección especificada:"

# ╔═╡ be11d027-1aef-4c3d-9f7e-01d70e5ca090
md"### 8.2.1  
$\vec{F} = x^2\hat{i} + 2x\hat{j} + z^2\hat{k}$. $C$: Elipse $4x^2 + y^2 = 4$ en el plano $xy$, en sentido antihorario visto desde arriba."

# ╔═╡ e393e493-1696-4d7a-b547-3c948baffb43
let
	F= [x^2,2x,z^2]
rotF=[
 diff(F[3], y) - diff(F[2], z),
 diff(F[1], z) - diff(F[3], x),
 diff(F[2], x) - diff(F[1], y)
]

end

# ╔═╡ d5b475c4-394e-4712-89b0-e2da77af557d
let
	F= [x^2,2x,z^2]
rotF=[
 diff(F[3], y) - diff(F[2], z),
 diff(F[1], z) - diff(F[3], x),
 diff(F[2], x) - diff(F[1], y)
]

g=[x,y,0]
∂gx∂gy=cross(diff.(g,x),diff.(g,y))
	dot(rotF,∂gx∂gy)
	
	integrate(dot(rotF,∂gx∂gy)*2r,(r,0,1),(t,0,2PI))
end

# ╔═╡ 4c86186f-1e41-443d-bbeb-a268556a7a02
md"### 8.2.2  
$\vec{F} = (x^2 + z^2)\hat{i} + (x^2 + y^2)\hat{j} + (x^2 + y^2)\hat{k}$. $C$: Cuadrado limitado por $x = \pm 1$ y $y = \pm 1$ en el plano $xy$, sentido antihorario visto desde arriba."

# ╔═╡ 52856b8f-4ebd-417e-a518-01841c98efec
let
	F= [x^2+z^2,x^2+y^2,x^2+y^2]
rotF=[
 diff(F[3], y) - diff(F[2], z),
 diff(F[1], z) - diff(F[3], x),
 diff(F[2], x) - diff(F[1], y)
]

g=[x,y,0]
∂gx∂gy=cross(diff.(g,x),diff.(g,y))
integrate(dot(rotF,∂gx∂gy),(x,-1,1),(y,-1,1))
end

# ╔═╡ 801b1d59-6c51-4256-b021-2ff27aa2289d
md"## 8.3 Integral del Curl sobre una Superficie"

# ╔═╡ c5074d09-baca-4b0f-8ad7-e794c6ffea99
md"### 8.3.1  
Si $\vec{n}$ es el vector unitario normal exterior de la superficie elíptica $S: 4x^2 + 9y^2 + 36z^2 = 36, \ z \ge 0$ y $\vec{F} = y\hat{i} + x^2\hat{j} + [x^2 + y^4]^{3/2} \sin(e^{\sqrt{xyz}})\hat{k}$.  
Encontrar el valor de $\iint_{S} (\nabla \times \vec{F}) \cdot \vec{n} \, d\sigma$.  
*Sugerencia: Usar la parametrización de la frontera $x = 3\cos(t), y = 2\sin(t), \ 0 \le t \le 2\pi$.*"

# ╔═╡ bb272aad-1b16-4de0-8f72-eea49becf49e
let
	F= [y,x^2,(x^2 +y^4)^(3/2) *sin(exp(sqrt(x*y*z)))]
rotF=[
 diff(F[3], y) - diff(F[2], z),
 diff(F[1], z) - diff(F[3], x),
 diff(F[2], x) - diff(F[1], y)
]

g=[x,y,0]
∂gx∂gy=cross(diff.(g,x),diff.(g,y))
dot(rotF,∂gx∂gy)
end

# ╔═╡ 5f1df999-f557-4be3-b437-76dd4e542c6f
let
	F= [y,x^2,(x^2 +y^4)^(3/2) *sin(exp(sqrt(x*y*z)))]
rotF=[
 diff(F[3], y) - diff(F[2], z),
 diff(F[1], z) - diff(F[3], x),
 diff(F[2], x) - diff(F[1], y)
]

g=[x,y,0]
∂gx∂gy=cross(diff.(g,x),diff.(g,y))
dot(rotF,∂gx∂gy)
	integrate(2*3*cos(t)-1*6r,(r,0,1),(t,0,2PI))
end

# ╔═╡ 9eee0fd8-d9e7-4022-aa4d-3e9994e87c2f
md"### 8.3.2  
Si $S$ es el cilindro $x^2 + y^2 = a^2, \ 0 \le z \le h$ junto con su tapa $x^2 + y^2 \le a^2, \ z = h$.  
Considerar $\vec{F} = -y\hat{i} + x\hat{j} + x^2\hat{k}$. Emplear el Teorema de Stokes para encontrar el flujo de $\nabla \times \vec{F}$ hacia afuera a través de $S$."

# ╔═╡ 3e7dad68-ad56-4d4e-9106-5138172d6fe5
let
	F= [-y,x,x^2]
rotF=[
 diff(F[3], y) - diff(F[2], z),
 diff(F[1], z) - diff(F[3], x),
 diff(F[2], x) - diff(F[1], y)
]

g=[r*cos(t),r*sin(t),0]
∂gx∂gy=cross(diff.(g,t),diff.(g,r))
dot(rotF,∂gx∂gy)
	integrate(dot(rotF,∂gx∂gy),(r,0,a),(t,0,2PI))
end

# ╔═╡ a6ba82a7-ae0d-4bb7-a97b-3eb3c64fdff8
md"## 8.4 Flujo del Curl a través de Superficies Parametrizadas
Calcular el flujo del curl de $\vec{F}$ en la dirección del vector normal unitario hacia afuera."

# ╔═╡ 87d7c868-1bf0-417c-9cce-db2d2eb27904
md"### 8.4.1  
$\vec{F} = 2z\hat{i} + 3x\hat{j} + 5y\hat{k}$.  
$S: \vec{r}(r, \theta) = r\cos\theta\hat{i} + r\sin\theta\hat{j} + (4 - r^2)\hat{k}, \ 0 \le r \le 2, \ 0 \le \theta \le 2\pi$."

# ╔═╡ c5ff98bd-65b4-48c8-a263-acaa58485adc
let
	F= [2z,3x,5y]
rotF=[
 diff(F[3], y) - diff(F[2], z),
 diff(F[1], z) - diff(F[3], x),
 diff(F[2], x) - diff(F[1], y)
]

g=[r*cos(t),r*sin(t),4-r^2]
∂gx∂gy=cross(diff.(g,t),diff.(g,r))
dot(rotF,∂gx∂gy)
	integrate(dot(rotF,∂gx∂gy),(r,0,2),(t,0,2PI))
end

# ╔═╡ 95664fd8-be27-4476-8fbe-4c162685a70f
md"### 8.4.2  
$\vec{F} = x^2y\hat{i} + 2y^3z\hat{j} + 3z\hat{k}$.  
$S: \vec{r}(r, \theta) = r\cos\theta\hat{i} + r\sin\theta\hat{j} + r\hat{k}, \ 0 \le r \le 1, \ 0 \le \theta \le 2\pi$."

# ╔═╡ 09753a66-d474-4f0a-b295-7dd73bc320e0
let
	F= [x^2*y,2y^3*z,3*z]
rotF=[
 diff(F[3], y) - diff(F[2], z),
 diff(F[1], z) - diff(F[3], x),
 diff(F[2], x) - diff(F[1], y)
]

g=[r*cos(t),r*sin(t),r]
	
∂gx∂gy=simplify.(cross(diff.(g,r),diff.(g,t)))

	integrate(dot([-2(r*sin(t))^3,0,-(r*cos(t))^2],∂gx∂gy),(r,0,1),(t,0,2PI))
end

# ╔═╡ 881322ce-b6d8-4637-9238-d39c98d160e4
md"### 8.4.3  
$\vec{F} = 3y\hat{i} + (5 - 2x)\hat{j} + (z^2 - 2)\hat{k}$.  
$S: \vec{r}(\phi, \theta) = \sqrt{3}\sin\phi\cos\theta\hat{i} + \sqrt{3}\sin\phi\sin\theta\hat{j} + \sqrt{3}\cos\phi\hat{k}$"

# ╔═╡ 807850a4-7c3f-475f-aee6-2dfc53ba378c
let
	F= [2y,(5-2x),z^2-2]
rotF=[
 diff(F[3], y) - diff(F[2], z),
 diff(F[1], z) - diff(F[3], x),
 diff(F[2], x) - diff(F[1], y)
]

g=[sqrt(3)*sin(p)*cos(t),sqrt(3)*sin(p)*sin(t),sqrt(3)*cos(p)]
∂gx∂gy=cross(diff.(g,p),diff.(g,t))
dot(rotF,∂gx∂gy)
	integrate(dot(rotF,∂gx∂gy),(p,0,sqrt(3)),(t,0,2PI))
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
# ╠═c6f01dbc-9a83-492e-9389-16776fb76c3d
# ╠═7c210c25-07b9-400b-b75e-bd4c9c8d5e68
# ╠═2be2dcd5-987b-4677-86a6-782f9ce6ff02
# ╟─d34d4f7a-1832-4b18-a6d0-01c41063146a
# ╟─eedce593-488f-4707-870c-ebc33ea96641
# ╟─001176ce-cff7-49ff-b1d7-47a40885faf8
# ╠═ea8109f3-7931-42be-bf42-cc994fdcf62c
# ╟─f1ecd628-a2a8-4daa-bd84-e6435f320935
# ╠═d04316c5-964a-478e-b9b4-6e8f64cf8115
# ╟─42c94328-3d8b-486a-9fed-8378171904dc
# ╠═4d8c93ba-a26c-488f-a221-57e831a02daf
# ╟─a23a72c8-5862-45e4-8512-19545b21bd48
# ╟─be11d027-1aef-4c3d-9f7e-01d70e5ca090
# ╠═e393e493-1696-4d7a-b547-3c948baffb43
# ╠═d5b475c4-394e-4712-89b0-e2da77af557d
# ╟─4c86186f-1e41-443d-bbeb-a268556a7a02
# ╠═52856b8f-4ebd-417e-a518-01841c98efec
# ╟─801b1d59-6c51-4256-b021-2ff27aa2289d
# ╟─c5074d09-baca-4b0f-8ad7-e794c6ffea99
# ╠═bb272aad-1b16-4de0-8f72-eea49becf49e
# ╠═5f1df999-f557-4be3-b437-76dd4e542c6f
# ╟─9eee0fd8-d9e7-4022-aa4d-3e9994e87c2f
# ╠═3e7dad68-ad56-4d4e-9106-5138172d6fe5
# ╟─a6ba82a7-ae0d-4bb7-a97b-3eb3c64fdff8
# ╟─87d7c868-1bf0-417c-9cce-db2d2eb27904
# ╠═c5ff98bd-65b4-48c8-a263-acaa58485adc
# ╟─95664fd8-be27-4476-8fbe-4c162685a70f
# ╠═09753a66-d474-4f0a-b295-7dd73bc320e0
# ╟─881322ce-b6d8-4637-9238-d39c98d160e4
# ╠═807850a4-7c3f-475f-aee6-2dfc53ba378c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
