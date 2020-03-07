var documenterSearchIndex = {"docs":
[{"location":"pages/api/#Functions-exported-from-ModifiedHankelFunctionsOfOrderOneThird:-1","page":"API","title":"Functions exported from ModifiedHankelFunctionsOfOrderOneThird:","text":"","category":"section"},{"location":"pages/api/#","page":"API","title":"API","text":"Modules = [ModifiedHankelFunctionsOfOrderOneThird]\nPrivate = false","category":"page"},{"location":"pages/api/#ModifiedHankelFunctionsOfOrderOneThird.modifiedhankel-Tuple{Any}","page":"API","title":"ModifiedHankelFunctionsOfOrderOneThird.modifiedhankel","text":"modifiedhankel(z)\n\nReturn h₁, h₂, h₁, and h₂, the first and second modified Hankel functions of order 1/3 and their derivatives.\n\nThese functions solve Stokes' equation d²udz² + zu = 0 as a power series in z for abs2(z) < 36 and an approximate asymptotic expansion otherwise. The asymptotic solution is necessary because the z³ⁱ in the power series blows up as i  .\n\nFor more information, see [SCL1945].\n\nSee also: powerseries, asymptotic\n\nExamples\n\njulia> h1, h2, h1prime, h2prime = modifiedhankel(complex(2.687, -0.648));\n\nReferences\n\n[SCL1945]: The Staff of the Computation Library (1945), Tables of the modified Hankel function of order one-third and of their derivatives. Cambridge, MA: Harvard University Press.\n\n\n\n\n\n","category":"method"},{"location":"pages/api/#Private-functions:-1","page":"API","title":"Private functions:","text":"","category":"section"},{"location":"pages/api/#","page":"API","title":"API","text":"Modules = [ModifiedHankelFunctionsOfOrderOneThird]\nPublic = false","category":"page"},{"location":"pages/api/#ModifiedHankelFunctionsOfOrderOneThird.asymptotic-Union{Tuple{T}, Tuple{T}} where T","page":"API","title":"ModifiedHankelFunctionsOfOrderOneThird.asymptotic","text":"asymptotic(z)\n\nReturn h₁, h₂, h₁, and h₂, the first and second modified Hankel functions of order 1/3 and their derivatives using an asymptotic expansion, useful for abs2(z) > 36.\n\nFor more information, see [SCL1945].\n\nSee also: modifiedhankel, powerseries\n\nReferences\n\n[SCL1945]: The Staff of the Computation Library (1945), Tables of the modified Hankel function of order one-third and of their derivatives. Cambridge, MA: Harvard University Press.\n\n\n\n\n\n","category":"method"},{"location":"pages/api/#ModifiedHankelFunctionsOfOrderOneThird.powerseries-Tuple{Any}","page":"API","title":"ModifiedHankelFunctionsOfOrderOneThird.powerseries","text":"powerseries(z)\n\nReturn h₁, h₂, h₁, and h₂, the first and second modified Hankel functions of order 1/3 and their derivatives using a power series with 30 terms.\n\nFor more information, see [SCL1945].\n\nSee also: modifiedhankel, asymptotic\n\nReferences\n\n[SCL1945]: The Staff of the Computation Library (1945), Tables of the modified Hankel function of order one-third and of their derivatives. Cambridge, MA: Harvard University Press.\n\n\n\n\n\n","category":"method"},{"location":"#Modified-Hankel-Functions-of-Order-One-Third-1","page":"Home","title":"Modified Hankel Functions of Order One Third","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Calculate modified Hankel functions of order one-third and their derivatives.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"These special functions are solutions to Stokes' differential equation:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"fracmathrmd^2umathrmdz^2 + zu = 0","category":"page"},{"location":"#","page":"Home","title":"Home","text":"The modified Hankel functions are linearly related to Airy functions but are not as well known. One application of the modified Hankel functions of order one-third is as a solution to the wave equation for electromagnetic fields between the boundaries of the earth-ionosphere waveguide. It may also occur in other instances of diffraction and refraction of waves.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"The direct computation of solutions to Stokes' equation is preferred to using Airy or Hankel functions. Unlike Bessel's equation, Stokes' equation has no singularity in the finite complex plane and its solutions are single-valued [SCL1945].","category":"page"},{"location":"#Usage-1","page":"Home","title":"Usage","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"using ModifiedHankelFunctionsOfOrderOneThird\n\nh1, h2, h1prime, h2prime = modifiedhankel(z)","category":"page"},{"location":"#Solutions-1","page":"Home","title":"Solutions","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Two solution approaches are used, as presented in [SCL1945]. If abs2(z) < 36, a power series solution is used. Otherwise, an asymptotic expansion is performed because of floating point limits in the power series.","category":"page"},{"location":"#Power-Series-Solution-1","page":"Home","title":"Power Series Solution","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"h₁, h₂, h₁, and h₂ are computed from auxiliary functions","category":"page"},{"location":"#","page":"Home","title":"Home","text":"h_1(z) = g + fracisqrt33(g - 2f) qquad h_2(z) = g - fracisqrt33(g - 2f)","category":"page"},{"location":"#","page":"Home","title":"Home","text":"where","category":"page"},{"location":"#","page":"Home","title":"Home","text":"f(z) = a_0 + a_1z^3 + a_2z^6 + cdots + a_mz^3m + cdots","category":"page"},{"location":"#","page":"Home","title":"Home","text":"g(z) = z(b_0 + b_1z^3 + b_2z^6 + cdots + b_mz^3m + cdots)","category":"page"},{"location":"#","page":"Home","title":"Home","text":"f(z) = -z^2(c_0 + c_1z^3 + c_2z^6 + cdots + c_mz^3m + cdots)","category":"page"},{"location":"#","page":"Home","title":"Home","text":"g(z) = d_0 + d_1z^3 + d_2z^6 + cdots + d_mz^3m + cdots","category":"page"},{"location":"#","page":"Home","title":"Home","text":"where","category":"page"},{"location":"#","page":"Home","title":"Home","text":"a_0 = frac2^13Gammaleft(frac23right) qquad a_m = -fraca_m-1(3m-1)3m","category":"page"},{"location":"#","page":"Home","title":"Home","text":"b_0 = frac2^133^23Gammaleft(frac43right) qquad b_m = -fracb_m-13m(3m+1)","category":"page"},{"location":"#","page":"Home","title":"Home","text":"c_0 = fraca_02 qquad c_m = - fracc_m-13m(3m+2)","category":"page"},{"location":"#","page":"Home","title":"Home","text":"d_0 = b_0 qquad d_m = - fracd_m-1(3m-2)3m","category":"page"},{"location":"#Asymptotic-Solution-1","page":"Home","title":"Asymptotic Solution","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"For h₁ on -2π3  arg z  4π3:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"h₁(z)  α z^-14 exp(23 i z^32 - 5πi12) left( 1 + sum_m=1 (-i)^m C_m z^-3m2 right)","category":"page"},{"location":"#","page":"Home","title":"Home","text":"h₁(z)  α i z^14 exp(23 i z^32 - 5πi12) left( 1 + sum_m=1 (-i)^m C_m z^-3m2 right) \n    - α4 z^-54 exp(23 i z^32 - 5πi12) left( 1 + sum_m=1 (-i)^m C_m z^-3m2 right) \n    - 32 α z^-14 exp(23 i z^32 - 5πi12) left( sum_m=1 (-i)^m m C_m z^-3m2 - 1 right)","category":"page"},{"location":"#","page":"Home","title":"Home","text":"where","category":"page"},{"location":"#","page":"Home","title":"Home","text":"C_m = frac(9-4)(81-4)cdots (92m-1^2-4)2^4m3^m m","category":"page"},{"location":"#","page":"Home","title":"Home","text":"and","category":"page"},{"location":"#","page":"Home","title":"Home","text":"h₁(z)  h₁(z) + α z^-14 exp(-23 i z^32 - 11πi12) left( 1 + sum_m_1 (i)^m C_m z^-3m2 right)","category":"page"},{"location":"#","page":"Home","title":"Home","text":"h₁(z)  h₁(z) - α z^14 exp(-23 i z^32 - 11πi12) left( 1 + sum_m_1 (i)^m C_m z^-3m2 right) \n    - α4 z^-54 exp(-23 i z^32 - 11πi12) left( 1 + sum_m_1 (i)^m C_m z^-3m2 right) \n    - 32 α z^-14 exp(-23 i z^32 - 11πi12) left( sum_m=1 i^m m C_m z^-3m2 - 1 right)","category":"page"},{"location":"#","page":"Home","title":"Home","text":"for -4π3  arg z  0.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"And for h₂ on -4π3  arg z  2π3:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"h₂(z)  α z^-14 exp(-23 i z^32 + 5πi12) left( 1 + sum_m=1 (i)^m C_m z^-3m2 right)","category":"page"},{"location":"#","page":"Home","title":"Home","text":"h₂(z)  -α z^14 exp(-23 i z^32 + 5πi12) left( 1 + sum_m=1 (i)^m C_m z^-3m2 right) \n    - α4 z^-54 exp(-23 i z^32 + 5πi12) left( 1 + sum_m=1 (i)^m C_m z^-3m2 right) \n    - 32 α z^-14 exp(-23 i z^32 + 5πi12) left( i^m m C_m z^-3m2 - 1 right)","category":"page"},{"location":"#","page":"Home","title":"Home","text":"and","category":"page"},{"location":"#","page":"Home","title":"Home","text":"h₂(z)  h₂(z) + α z^-14 exp(23 i z^32 + 11πi12) left( 1 + sum_m=1 (-i)^m C_m z^-3m2 right)","category":"page"},{"location":"#","page":"Home","title":"Home","text":"h₂(z)  h₂(z) + α z^14 exp(23 i z^32 + 11πi12) left( 1 + sum_m=1 (-i)^m C_m z^-3m2 right) \n    - α4 z^-54 exp(23 i z^32 + 11πi12) left( 1 + sum_m=1 (-i)^m C_m z^-3m2 right) \n    - 32 α z^-14 exp(23 i z^32 + 11πi12) left( sum_m=1 (-i)^m m C_m z^-3m2 - 1 right)","category":"page"},{"location":"#","page":"Home","title":"Home","text":"for 0  arg z  4π3.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"The coefficient alpha = 2^13 3^16 pi^-12.","category":"page"},{"location":"#References-1","page":"Home","title":"References","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"[SCL1945]: The Staff of the Computation Library (1945), Tables of the modified Hankel function of order one-third and of their derivatives. Cambridge, MA: Harvard University Press.","category":"page"}]
}
