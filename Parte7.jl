### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ ae48d93d-79a1-4453-b582-902245c6f310
using SymPy

# ╔═╡ d65b6161-dd72-4b71-ac50-ff35a5ad499f
using LinearAlgebra

# ╔═╡ da17e671-5e2e-4709-833e-34770d68e49d
md"# 7 Integrales de Superficie"

# ╔═╡ 780bc55d-8b86-4984-915a-771a4f26ef0a
@syms x y z t θ::real r ϕ::real u v c a::real b c p

# ╔═╡ 686b7daa-6758-45a0-b7b4-efe8d53f4702
md"## 7.1 Integrales de Superficie de Funciones Escalares
Encontrar la integral de superficie $\iint_{S} G(x, y, z) \, d\sigma$ para las siguientes funciones y superficies:"

# ╔═╡ c8627ae6-7097-41ae-ab8b-bd267d451411
md"### 7.1.1  
$G(x, y, z) = x$ sobre el cilindro parabólico $y = x^2, \quad 0 \le x \le 2, \quad 0 \le z \le 3$"

# ╔═╡ 0d417090-45f8-4dfe-bc82-9be949bcada6
let
	g=x
	f=x^2
    fx=diff(f,x)
 fy=diff(f,y)
	
	integrate(g*sqrt((fx)^2 +(fy)^2 +1),(x,0,2),(z,0,3)).evalf()
	
end

# ╔═╡ ab1a2f8a-829e-4f8d-a626-8deeb2984ef5
md"### 7.1.2  
$F(x, y, z) = z$ sobre la porción del plano $x + y + z = 4$ que está sobre el cuadrado $0 \le x \le 1, \quad 0 \le y \le 1$ en el plano $xy$."

# ╔═╡ e507354e-07db-4349-ba18-b48da510927c
let
	g=4-x-y
	f=x+y+z
    fx=diff(f,x)
 fy=diff(f,y)
	

gF=norm([diff(f,x),diff(f,y),diff(f,z)])
	gF_k=norm(dot([diff(f,x),diff(f,y),diff(f,z)],[0,0,1]))

	integrate(g*gF/gF_k,(x,0,1),(y,0,1)).evalf()
	
	
end

# ╔═╡ d5f7410f-3b13-4380-bec7-442d49420962
let
	g=4-x-y
	f=4-x-y
    fx=diff(f,x)
 fy=diff(f,y)
	
	integrate(g*sqrt((fx)^2 +(fy)^2 +1),(x,0,1),(y,0,1)).evalf()
	
end

# ╔═╡ 455bd615-ee74-436c-9dca-9bcdc119f13e
md"### 7.1.3  
$H(x, y, z) = x^2 \sqrt{5 - 4z}$ sobre el domo parabólico $z = 1 - x^2 - y^2, \quad z \ge 0$"

# ╔═╡ 501fa9fa-258d-4759-a62e-fcb339fdca8f
let
	g=x^2*sqrt(5-4*z)
	f=1-x^2-y^2
    fx=diff(f,x)
 fy=diff(f,y)
	
	g*sqrt((fx)^2 +(fy)^2 +1)
	
end

# ╔═╡ 7c5a3866-4e7b-4f04-ac5a-14fc82244e7e
let
	f=(r*cos(θ))^2*sqrt(5-4(1-r^2))*sqrt(4*r^2 + 1)

	integrate(f*r,(r,0,1),(θ,0,2PI))
end

# ╔═╡ 0005958d-d5bc-4451-afbf-203d10876750
md"## 7.2 Superficie de un Sólido Rectangular
Integrar $G(x, y, z) = xyz$ sobre la superficie del sólido rectangular recortado en el primer octante por los planos $x = a, \quad y = b, \quad z = c$."

# ╔═╡ 5d438ddc-1eec-468f-86d5-fe0e83e0364e
let
	g=x*y*c
	f=c
    fx=diff(f,x)
 fy=diff(f,y)
	
	integrate(g*sqrt((fx)^2 +(fy)^2 +1),(x,0,a),(y,0,b)).evalf()
	
end

# ╔═╡ 708f323d-e692-4d8d-9c75-229cd109085a
let
	g=x*b*z
	f=b
    fx=diff(f,x)
 fy=diff(f,z)
	
	integrate(g*sqrt((fx)^2 +(fy)^2 +1),(x,0,a),(z,0,c)).evalf()
	
end

# ╔═╡ e18805bb-4046-4e0a-9dcc-daa6a084f48f
let
	g=a*y*z
	f=a
    fx=diff(f,x)
 fy=diff(f,z)
	
	integrate(g*sqrt((fx)^2 +(fy)^2 +1),(y,0,b),(z,0,c)).evalf()
	
end

# ╔═╡ b97c3416-396a-47ac-ad2b-828e58b4e3b1
let
	t=0
	g=a*y*z
	f=a
    fx=diff(f,x)
 fy=diff(f,z)
	
	t=integrate(g*sqrt((fx)^2 +(fy)^2 +1),(y,0,b),(z,0,c)).evalf()+t
	
	g=x*y*c
	f=c
    fx=diff(f,x)
 fy=diff(f,y)
	
	t=integrate(g*sqrt((fx)^2 +(fy)^2 +1),(x,0,a),(y,0,b)).evalf()+t

		g=a*y*z
	f=a
    fx=diff(f,x)
 fy=diff(f,z)
	
	t=integrate(g*sqrt((fx)^2 +(fy)^2 +1),(y,0,b),(z,0,c)).evalf()+t
end

# ╔═╡ da75e0eb-3c57-4ab7-af6c-ee0f3a84b35d
md"## 7.3 Superficie Triangular
Integrar $G(x, y, z) = xyz$ sobre la superficie triangular con vértices $(1, 0, 0), \quad (0, 2, 0), \quad$ y $(0, 1, 1)$."

# ╔═╡ 3ac75703-6408-4b56-9f9f-87a3393fa15a
let
    x = 1 - u - v
    y = 2*u + v
    z = v
    
    g = x * y * z
    
     xu = diff(x, u)
    yu = diff(y, u)
      zu = diff(z, u)
    
     
	
	xv = diff(x, v)
    yv = diff(y, v)
     zv = diff(z, v)

    cx = yu*zv - zu*yv
    cy = zu*xv - xu*zv
    cz = xu*yv - yu*xv
    mag = sqrt(cx^2 + cy^2 + cz^2)
    
  integrate(g * mag, (v, 0, 1 - u), (u, 0, 1))
   
end

# ╔═╡ 75b9d303-25cc-4577-a245-9ed813e6b7ac
md"## 7.4 Flujo de Campos Vectoriales (Flux)
Usar una parametrización para encontrar el flujo $\iint_{S} \vec{F} \cdot \vec{n} \, d\sigma$ a través de la superficie en la dirección especificada."

# ╔═╡ 7231e2c3-d2a5-4ed7-a28f-226d28a33d69
md"### 7.4.1  
$\vec{F} = z^2 \hat{i} + x \hat{j} - 3z \hat{k}$ hacia afuera (normal alejándose del eje $x$) a través de la superficie formada por el corte de los planos $x = 0, \quad x = 1, \quad z = 0$ sobre el cilindro $z = 4 - y^2$."

# ╔═╡ 25a31499-3040-4b18-a6f1-b35a51916858
let
	r=[x,y,(4-y^2)]

	n=(cross(diff.(r,x),diff.(r,y)))
	
	integrate(dot([(4-y^2)^2,x,-3(4-y^2)],n),(x,0,1),(y,-2,2)).evalf()
end

# ╔═╡ 909a048d-1ec3-4666-a490-438c9fb3d5f7
md"### 7.4.2  
$\vec{F} = 2xy \hat{i} + 2yz \hat{j} + 2xz \hat{k}$ hacia arriba, a través de la porción del plano $x + y + z = 2a$ sobre el cuadrado $0 \le x \le a, \quad 0 \le y \le a$ en el plano $xy$."

# ╔═╡ ad3e2dd1-ac9f-4fa0-bf39-9fc40f1a7c42
let
r=[x,y,2*a-x-y]
n=(cross(diff.(r,x),diff.(r,y)))
Fn=dot([2*x*y,2*y*(2*a-x-y),2*x*(2*a-x-y)],n)

end

# ╔═╡ e541d74f-2231-496c-af78-ce53ff85fec1
	integrate(2(2*a-x-y)*x+2(2*a-x-y)*y+2*x*y,(x,0,a),(y,-0,a))

# ╔═╡ 67add63e-851a-4637-b27f-aad8060de426
md"### 7.4.3  
$\vec{F} = xy \hat{i} - z \hat{k}$ hacia afuera (normal alejándose del eje $z$) a través del cono $z = \sqrt{x^2 + y^2}, \quad 0 \le z \le 1$."

# ╔═╡ 22997efd-5c8f-4053-be77-c163f24d571d
let
R=[r*cos(θ),r*sin(θ),r]
n=(-1)(cross(diff.(R,r),diff.(R,θ)))
Fn=dot([r^2*sin(θ)*cos(θ),0,-r],n)

end

# ╔═╡ 87bc93a1-a274-494f-aa7f-26cf7c8b19c8
integrate(r^3*cos(θ)^2*sin(θ) + r^2, (r, 0, 1), (θ, 0, 2*PI))

# ╔═╡ 87217b67-d3cd-44ef-bd30-220cddc1c851
md"### 7.4.4  
$\vec{F} = -x \hat{i} - y \hat{j} + z^2 \hat{k}$ hacia afuera (normal alejándose del eje $z$) a través de la porción de cono $z = \sqrt{x^2 + y^2}$ entre los planos $z = 1$ y $z = 2$."

# ╔═╡ f24ea474-887f-41e5-9f48-6403a2fcc89c
let
R=[r*cos(θ),r*sin(θ),r]
n=(cross(diff.(R,θ),diff.(R,r)))
Fn=dot([-r*cos(θ),-r*sin(θ),r^2],n)
F_n_do=-r^2*sin(θ)^2-r^2*cos(θ)^2-r^3

end

# ╔═╡ 63b827b4-1508-4eb7-939b-9d82a702a93b
let
R=[r*cos(θ),r*sin(θ),r]
n=(cross(diff.(R,θ),diff.(R,r)))
Fn=dot([-r*cos(θ),-r*sin(θ),r^2],n)
F_n_do=-r^2*sin(θ)^2-r^2*cos(θ)^2-r^3
integrate((F_n_do), (r, 1, 2), (θ, 0, 2*PI))
end

# ╔═╡ 162589fc-5477-4f91-9790-2596a8184b30
md"## 7.5 Flujo sobre una Esfera
Encontrar el Flux empleando la integral de superficie del campo vectorial $\vec{F}$ sobre la porción de esfera $x^2 + y^2 + z^2 = a^2$ en el primer octante, en la dirección alejándose del origen."

# ╔═╡ 4a51f12e-2037-4782-b89a-0d1fd6548219
md"### 7.5.1  
$\vec{F}(x, y, z) = z\hat{k}$"

# ╔═╡ d662d20e-c0aa-4a27-b655-1c5a3562a3ab
integrate(p^2*sin(θ),(p,0,a),(θ,0,PI/2),(ϕ,0,PI/2))

# ╔═╡ f7bed39d-0b3c-4079-bdcf-94a58630d7d9
md"### 7.5.2  
$\vec{F}(x, y, z) = y\hat{i} - x\hat{j} + \hat{k}$"

# ╔═╡ 1258a09c-1c08-4b25-a4ac-e705ab939e0b
let
R=[a*cos(θ)*sin(ϕ),a*sin(ϕ)sin(θ),a*cos(ϕ)]
rx=diff.(R,ϕ)
ry=diff.(R,θ)
	
	f=[a*sin(ϕ)sin(θ),-a*cos(θ)*sin(ϕ),1]
	
	integrate(dot(f,cross(rx,ry)),(ϕ,0,PI/2),(θ,0,PI/2))
end

# ╔═╡ cf24e461-f232-45eb-9db0-03175701addd
md"### 7.5.3  
$\vec{F}(x, y, z) = x\hat{i} + y\hat{j} + z\hat{k}$"

# ╔═╡ aaf9e8d2-5563-4a84-a961-f161ce1e921d
let
R=[a*cos(θ)*sin(ϕ),a*sin(ϕ)sin(θ),a*cos(ϕ)]
rx=diff.(R,ϕ)
ry=diff.(R,θ)
	
	f=[a*cos(θ)*sin(ϕ),a*sin(ϕ)sin(θ),a*cos(ϕ)]
	
	integrate(dot(f,cross(rx,ry)),(ϕ,0,PI/2),(θ,0,PI/2))
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
# ╠═ae48d93d-79a1-4453-b582-902245c6f310
# ╠═d65b6161-dd72-4b71-ac50-ff35a5ad499f
# ╟─da17e671-5e2e-4709-833e-34770d68e49d
# ╠═780bc55d-8b86-4984-915a-771a4f26ef0a
# ╟─686b7daa-6758-45a0-b7b4-efe8d53f4702
# ╟─c8627ae6-7097-41ae-ab8b-bd267d451411
# ╠═0d417090-45f8-4dfe-bc82-9be949bcada6
# ╟─ab1a2f8a-829e-4f8d-a626-8deeb2984ef5
# ╠═e507354e-07db-4349-ba18-b48da510927c
# ╠═d5f7410f-3b13-4380-bec7-442d49420962
# ╟─455bd615-ee74-436c-9dca-9bcdc119f13e
# ╠═501fa9fa-258d-4759-a62e-fcb339fdca8f
# ╠═7c5a3866-4e7b-4f04-ac5a-14fc82244e7e
# ╟─0005958d-d5bc-4451-afbf-203d10876750
# ╠═5d438ddc-1eec-468f-86d5-fe0e83e0364e
# ╠═708f323d-e692-4d8d-9c75-229cd109085a
# ╠═e18805bb-4046-4e0a-9dcc-daa6a084f48f
# ╠═b97c3416-396a-47ac-ad2b-828e58b4e3b1
# ╟─da75e0eb-3c57-4ab7-af6c-ee0f3a84b35d
# ╠═3ac75703-6408-4b56-9f9f-87a3393fa15a
# ╟─75b9d303-25cc-4577-a245-9ed813e6b7ac
# ╟─7231e2c3-d2a5-4ed7-a28f-226d28a33d69
# ╠═25a31499-3040-4b18-a6f1-b35a51916858
# ╟─909a048d-1ec3-4666-a490-438c9fb3d5f7
# ╠═ad3e2dd1-ac9f-4fa0-bf39-9fc40f1a7c42
# ╠═e541d74f-2231-496c-af78-ce53ff85fec1
# ╟─67add63e-851a-4637-b27f-aad8060de426
# ╠═22997efd-5c8f-4053-be77-c163f24d571d
# ╠═87bc93a1-a274-494f-aa7f-26cf7c8b19c8
# ╟─87217b67-d3cd-44ef-bd30-220cddc1c851
# ╠═f24ea474-887f-41e5-9f48-6403a2fcc89c
# ╠═63b827b4-1508-4eb7-939b-9d82a702a93b
# ╟─162589fc-5477-4f91-9790-2596a8184b30
# ╟─4a51f12e-2037-4782-b89a-0d1fd6548219
# ╠═d662d20e-c0aa-4a27-b655-1c5a3562a3ab
# ╟─f7bed39d-0b3c-4079-bdcf-94a58630d7d9
# ╠═1258a09c-1c08-4b25-a4ac-e705ab939e0b
# ╟─cf24e461-f232-45eb-9db0-03175701addd
# ╠═aaf9e8d2-5563-4a84-a961-f161ce1e921d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
