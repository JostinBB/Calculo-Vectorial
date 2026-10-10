### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 275168eb-fc03-4475-b21f-1a895421816d
using SymPy

# ╔═╡ 2cec8a64-edde-4145-afb6-f8695a040bb6
using LinearAlgebra

# ╔═╡ 0a7409db-75d9-4f1b-a2a6-96443fe027f0
md"# 5 Teorema de Green en el Plano"

# ╔═╡ dadca21a-aa94-4b6b-9c81-875eb50559a8
@syms x y z t a r theta

# ╔═╡ d9cf6e7d-5a87-4d44-a8d7-eb5335831a93
md"## 5.1 Rotacional en el Plano
Encontrar la componente $k$ de $\text{rot} \, \vec{F}$ (curl $\vec{F}$) para los siguientes campos vectoriales:"

# ╔═╡ 5e4afa7e-9cf9-4f81-8de5-c11555bc7552
md"### 5.1.1  
$\vec{F} = (x^2 + y)\hat{i} + 2xy\hat{j}$"

# ╔═╡ 8724cef0-ec04-47f8-b05c-9c50d7723ba2
let
    P = x^2 + y
    Q = 2*x*y
    
    rotK = diff(Q, x) - diff(P, y)
end

# ╔═╡ 379c85ec-849f-43aa-a73e-f82806b17a74
md"### 5.1.2  
$\vec{F} = (xe^y)\hat{i} + (ye^x)\hat{j}$"

# ╔═╡ dd5e029a-01d9-4dd8-b285-f43666812004
let
    P = x*exp(y)
    Q = y*exp(x)
    
    rotK = diff(Q, x) - diff(P, y)
end

# ╔═╡ 9072856a-c477-4f8d-8611-621a8b0266ae
md"### 5.1.3  
$\vec{F} = (y \sin x)\hat{i} + (x \sin y)\hat{j}$"

# ╔═╡ 366c70d6-6e65-4548-8d5f-e5621ef54f71
let
    P = y*sin(x)
    Q = x*sin(y)
    
    rotK = diff(Q, x) - diff(P, y)
end

# ╔═╡ 12faf561-3786-45e4-99d9-fc637763b460
md"### 5.1.4  
$\vec{F} = \left( \frac{x}{y} \right)\hat{i} - \left( \frac{y}{x} \right)\hat{j}$"

# ╔═╡ cbd33db7-8a86-42bd-9373-d97b715b818e
let
    P = x/y
    Q = -y/x
    
    rotK = diff(Q, x) - diff(P, y)
end

# ╔═╡ b27c28bc-633a-417b-82fe-1fb997913e95
md"## 5.2 Verificación del Teorema de Green
Verificar el Teorema de Green evaluando tanto la integral de línea como la integral doble para obtener la **Circulación** y el **Flux**.  
*Dominio:* Disco $R: x^2 + y^2 \le a^2$.  
*Curva:* Círculo $C: \vec{r} = a\cos(t)\hat{i} + a\sin(t)\hat{j}, \quad 0 \le t \le 2\pi$."

# ╔═╡ 69ea4593-0bd0-4a31-a31e-5b3981541f92
md"circualcion"

# ╔═╡ 8bf46507-2b31-42e6-9c82-194e3e6c62b9
md"### 5.2.1  
$\vec{F} = -y\hat{i} + x\hat{j}$"

# ╔═╡ 2cf975ab-18eb-4950-a261-4e863f3bddcf
let
	integrate(dot([-a*sin(t),a*cos(t)],[-a*sin(t),a*cos(t)]),(t,0,2PI))
end

# ╔═╡ 40e2b5d0-0c28-4ceb-8352-57191c93bd56
# integrate(función, (variable1, min, max), (variable2, min, max))

# ╔═╡ 662cc181-790e-4e7f-8d63-a524db977ffd
let
	P=-y
	Q=x
	dQx=diff(Q,x)
	dPy=diff(P,y)
	dQx-dPy
end

# ╔═╡ 846072cd-03bc-44b0-b862-fa361fef0a9a
let
	P=-y
	Q=x
	dQx=diff(Q,x)
	dPy=diff(P,y)

	
	integrate(r*(dQx-dPy),(r,0,a),(t,0,2PI))
end

# ╔═╡ 313f2205-342f-4309-a55c-40accbdb6e58
md"flux"

# ╔═╡ f4418489-8550-47f3-9a4a-b593d58485d7
let
	P=-y
	Q=x
	dQx=diff(Q,y)
	dPy=diff(P,x)
	dQx+dPy
		integrate(r*(dQx+dPy),(r,0,a),(t,0,2PI))
end

# ╔═╡ 6f48eca4-55a9-4562-8e08-7d93418f07b6
md"### 5.2.2  
$\vec{F} = -x^2 y\hat{i} + xy^2 \hat{j}$"

# ╔═╡ 72335899-e740-47d7-97bb-e5e489631171
md"circualcion"

# ╔═╡ a8b65c9c-ac80-4900-85ba-7e952c44988e
let
	P=-x^2 *y
	Q=x*y^2
	dQx=diff(Q,x)
	dPy=diff(P,y)
	dQx-dPy
end

# ╔═╡ bff28649-7086-4773-8a04-621a2b93afb2
let
	
	integrate(r*(r^2),(r,0,a),(t,0,2PI))
end

# ╔═╡ efb8ce16-44b1-4739-848c-b06743864642
md"flux"

# ╔═╡ 77f5c64d-5616-4c0f-a501-aceff566feec
let
	P=-x^2 *y
	Q=x*y^2
	dQx=diff(Q,y)
	dPy=diff(P,x)
	dQx+dPy
end

# ╔═╡ d6e986b2-e7bc-4d81-92d8-cf1fc4c92e94
let
	P=-y
	Q=x
	dQx=diff(Q,y)
	dPy=diff(P,x)
	dQx+dPy
		integrate(r*(dQx+dPy),(r,0,a),(t,0,2PI))
end

# ╔═╡ 747432c2-4e8a-4972-89e1-9a6217c7e45c
md"## 5.3 Aplicación del Teorema de Green
Encontrar la circulación (sentido antihorario) y el flux hacia afuera para el campo $\vec{F}$ y la curva $C$:"

# ╔═╡ 23faac5b-b085-42f7-9156-062901b410f2
md"### 5.3.1  
$\vec{F} = (x - y)\hat{i} + (y - x)\hat{j}$. $C$: Cuadrado limitado por $x=0, x=1, y=0, y=1$."

# ╔═╡ 1b4247a9-9083-4b96-b3ab-ecaa912d08f2
let
	P=x-y
	Q=y-x
	dQx=diff(Q,x)
	dPy=diff(P,y)
	dQx-dPy
end

# ╔═╡ 32e191cb-706b-46b1-bbe8-9c75952442f7
md"como el rotacional es 0 la cirucalcion es 0"

# ╔═╡ 5df1d01b-4f4f-467f-ae66-2c0ab9cf568a
md"### 5.3.2  
$\vec{F} = (y^2 - x^2)\hat{i} + (x^2 + y^2)\hat{j}$. $C$: Triángulo limitado por $y=0, x=3, y=x$."

# ╔═╡ d29496c6-42ba-4e58-b45f-9d00cd429ce1
let
	P=y^2-x^2
	Q=y^2 + x^2
	dQx=diff(Q,x)
	dPy=diff(P,y)
	dQx-dPy

	integrate((dQx-dPy),(y,0,x),(x,0,3))
	
end

# ╔═╡ e40fbea8-754d-4892-ad43-224050485876
md"### 5.3.3  
$\vec{F} = (xy + y^2)\hat{i} + (x - y)\hat{j}$. $C$: Región encerrada por $y=x^2$ y $x=y^2$."

# ╔═╡ ac481cbe-adc1-45c1-add8-d9bae6b9e6e8
let
	P=x*y +y^2
	Q=x-y
	dQx=diff(Q,x)
	dPy=diff(P,y)
	dQx-dPy
	
end

# ╔═╡ 9617e0fe-faef-4b51-b89d-28da4a7a520b
let
	P=x*y +y^2
	Q=x-y
	dQx=diff(Q,x)
	dPy=diff(P,y)
	dQx-dPy
	integrate((dQx-dPy),(y,x^2,sqrt(x)),(x,0,1))
	
end

# ╔═╡ fad2c021-38a2-454d-97e5-3e9e120a504d
md"### 5.3.4  
$\vec{F} = (x + e^x \sin y)\hat{i} + (x + e^x \cos y)\hat{j}$. $C$: El bucle derecho de la curva polar $r^2 = \cos(2\theta)$."

# ╔═╡ 16e5f234-3d32-4972-9043-193e6c4d8859
let
	P=x+exp(x)*sin(y)
	Q=x+exp(x)*cos(y)
	dQx=diff(Q,x)
	dPy=diff(P,y)
	dQx-dPy
	2*integrate(r*(dQx-dPy),(r,0,sqrt(cos(2*t))),(t,0,PI/4))
end

# ╔═╡ 39b5ff37-8a51-4f78-8c46-3904165e3406
md"## 5.4 Flujo a través de un Cardioide
Encontrar el flux hacia afuera de $\vec{F}$ a través del cardioide $r = a(1 + \cos \theta), a > 0$.  
$\vec{F} = \left( 3xy - \frac{x}{1+y^2} \right)\hat{i} + (e^x + \tan^{-1} y)\hat{j}$"

# ╔═╡ 254e2b30-538e-477f-ad3a-dddec8e6a249
let
	P=3*x*y-(x/(1+y^2))
	Q=exp(x)+atan(y)
	dQx=diff(Q,y)
	dPy=diff(P,x)
	dQx+dPy
		
end

# ╔═╡ fcf23f3b-2abc-4d7d-899b-0922572d669c
let
	P=3*x*y-(x/(1+y^2))
	Q=exp(x)+atan(y)
	dQx=diff(Q,y)
	dPy=diff(P,x)
	dQx+dPy
			integrate(r*(3*r*sin(t)),(r,0,a*(1+cos(t))),(t,0,2*PI))
end

# ╔═╡ 8ef15a74-d987-41f6-a63f-2ef397b2d418
md"## 5.5 Trabajo realizado
Encontrar el trabajo para mover una partícula una vez en sentido antihorario.  
$\vec{F} = 2xy^3 \hat{i} + 4x^2 y^2 \hat{j}$. $C$: Borde de la región entre el eje $x$, $x=1$ e $y=x^3$."

# ╔═╡ 7a670c1a-a709-4e43-b5d7-6dce95f45876
let
	P=2*x*y^3
	Q=4*x^2*y^2
	dQx=diff(Q,x)
	dPy=diff(P,y)
	dQx-dPy

	
end

# ╔═╡ 3d41a981-1d0c-41ec-9260-fc01830f9702
let
	P=2*x*y^3
	Q=4*x^2*y^2
	dQx=diff(Q,x)
	dPy=diff(P,y)
	dQx-dPy
	integrate((dQx-dPy),(y,0,x^3),(x,0,1))
	
end

# ╔═╡ 0b096a3c-94f1-42ea-9ef6-d388d4259267
md"## 5.6 Integrales de Línea Cerrada
Evaluar usando el Teorema de Green:"

# ╔═╡ c94eee08-7035-4132-8899-84801323ea22
md"### 5.6.1  
$\oint_{C} y^2 \, dx + x^2 \, dy$. $C$: Triángulo definido por $x=0, x+y=1, y=0$."

# ╔═╡ 248d90d1-69dd-437b-99ab-efe1ef303d3f
let
	P=y^2
	Q=x^2
	dQx=diff(Q,x)
	dPy=diff(P,y)
	dQx-dPy
	integrate((dQx-dPy),(y,0,1-x),(x,0,1))
	
end

# ╔═╡ 49f58bd1-a31f-4b38-8b75-53fcd615e9a0
md"### 5.6.2  
$\oint_{C} (6y + x) \, dx + (y + 2x) \, dy$. $C$: Círculo $(x-2)^2 + (y-3)^2 = 4$."

# ╔═╡ 08c13a71-7985-4237-88d3-4852d7e188c8
let
	P=6*y+x
	Q=y+2*x
	dQx=diff(Q,x)
	dPy=diff(P,y)
	dQx-dPy

	integrate((dQx-dPy)*r,(r,0,2 ),(t,0,2*PI))
end

# ╔═╡ b53fdc10-5a76-4c6b-82f4-3b31739b5c14
md"## 5.7 Cálculo de Áreas con Teorema de Green
El área de una región $R$ encerrada por $C$ se define como:
$$\text{Área} = \frac{1}{2} \oint_{C} x \, dy - y \, dx$$

Encontrar el área encerrada por:"

# ╔═╡ cc3c8068-d0d8-4eae-95c0-8e19ebd5587d
md"### 5.7.1  
El círculo $\vec{r}(t) = a\cos(t)\hat{i} + a\sin(t)\hat{j}, \quad 0 \le t \le 2\pi$."

# ╔═╡ 26d8222c-de74-403b-8416-4a32c0b292a5
let
	P=-y
	Q=x
	dQx=diff(Q,x)
	dPy=diff(P,y)
	dQx-dPy
	
		1/2 *integrate((dQx-dPy)*r,(r,0,a ),(t,0,2*PI))
end

# ╔═╡ 4d39583e-89f5-433b-93fd-975ab4f41603
md"### 5.7.2  
El astroide $\vec{r}(t) = \cos^3(t)\hat{i} + \sin^3(t)\hat{j}, \quad 0 \le t \le 2\pi$."

# ╔═╡ 1061b530-c122-4a14-92aa-263455f301a8
let
	X=cos(t)^3
	Y=sin(t)^3
	dQx=diff(X,t)
	dPy=diff(Y,t)
	dQx-dPy
	1/2 *integrate(dot([-Y,X],[dQx,dPy]),(t,0,2*PI))
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
# ╠═275168eb-fc03-4475-b21f-1a895421816d
# ╠═2cec8a64-edde-4145-afb6-f8695a040bb6
# ╟─0a7409db-75d9-4f1b-a2a6-96443fe027f0
# ╠═dadca21a-aa94-4b6b-9c81-875eb50559a8
# ╟─d9cf6e7d-5a87-4d44-a8d7-eb5335831a93
# ╟─5e4afa7e-9cf9-4f81-8de5-c11555bc7552
# ╠═8724cef0-ec04-47f8-b05c-9c50d7723ba2
# ╟─379c85ec-849f-43aa-a73e-f82806b17a74
# ╠═dd5e029a-01d9-4dd8-b285-f43666812004
# ╟─9072856a-c477-4f8d-8611-621a8b0266ae
# ╠═366c70d6-6e65-4548-8d5f-e5621ef54f71
# ╟─12faf561-3786-45e4-99d9-fc637763b460
# ╠═cbd33db7-8a86-42bd-9373-d97b715b818e
# ╟─b27c28bc-633a-417b-82fe-1fb997913e95
# ╟─69ea4593-0bd0-4a31-a31e-5b3981541f92
# ╟─8bf46507-2b31-42e6-9c82-194e3e6c62b9
# ╠═2cf975ab-18eb-4950-a261-4e863f3bddcf
# ╟─40e2b5d0-0c28-4ceb-8352-57191c93bd56
# ╠═662cc181-790e-4e7f-8d63-a524db977ffd
# ╠═846072cd-03bc-44b0-b862-fa361fef0a9a
# ╟─313f2205-342f-4309-a55c-40accbdb6e58
# ╠═f4418489-8550-47f3-9a4a-b593d58485d7
# ╟─6f48eca4-55a9-4562-8e08-7d93418f07b6
# ╟─72335899-e740-47d7-97bb-e5e489631171
# ╠═a8b65c9c-ac80-4900-85ba-7e952c44988e
# ╠═bff28649-7086-4773-8a04-621a2b93afb2
# ╟─efb8ce16-44b1-4739-848c-b06743864642
# ╠═77f5c64d-5616-4c0f-a501-aceff566feec
# ╠═d6e986b2-e7bc-4d81-92d8-cf1fc4c92e94
# ╟─747432c2-4e8a-4972-89e1-9a6217c7e45c
# ╟─23faac5b-b085-42f7-9156-062901b410f2
# ╠═1b4247a9-9083-4b96-b3ab-ecaa912d08f2
# ╟─32e191cb-706b-46b1-bbe8-9c75952442f7
# ╟─5df1d01b-4f4f-467f-ae66-2c0ab9cf568a
# ╠═d29496c6-42ba-4e58-b45f-9d00cd429ce1
# ╟─e40fbea8-754d-4892-ad43-224050485876
# ╠═ac481cbe-adc1-45c1-add8-d9bae6b9e6e8
# ╠═9617e0fe-faef-4b51-b89d-28da4a7a520b
# ╟─fad2c021-38a2-454d-97e5-3e9e120a504d
# ╠═16e5f234-3d32-4972-9043-193e6c4d8859
# ╟─39b5ff37-8a51-4f78-8c46-3904165e3406
# ╠═254e2b30-538e-477f-ad3a-dddec8e6a249
# ╠═fcf23f3b-2abc-4d7d-899b-0922572d669c
# ╟─8ef15a74-d987-41f6-a63f-2ef397b2d418
# ╠═7a670c1a-a709-4e43-b5d7-6dce95f45876
# ╠═3d41a981-1d0c-41ec-9260-fc01830f9702
# ╟─0b096a3c-94f1-42ea-9ef6-d388d4259267
# ╟─c94eee08-7035-4132-8899-84801323ea22
# ╠═248d90d1-69dd-437b-99ab-efe1ef303d3f
# ╟─49f58bd1-a31f-4b38-8b75-53fcd615e9a0
# ╠═08c13a71-7985-4237-88d3-4852d7e188c8
# ╟─b53fdc10-5a76-4c6b-82f4-3b31739b5c14
# ╟─cc3c8068-d0d8-4eae-95c0-8e19ebd5587d
# ╠═26d8222c-de74-403b-8416-4a32c0b292a5
# ╟─4d39583e-89f5-433b-93fd-975ab4f41603
# ╠═1061b530-c122-4a14-92aa-263455f301a8
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
