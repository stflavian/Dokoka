using Dokoka

file = open("mock-build.xyz", "w")
write(file, 
"""12
Benzene
H      1.2194     -0.1652      2.1600
C      0.6825     -0.0924      1.2087
C     -0.7075     -0.0352      1.1973
H     -1.2644     -0.0630      2.1393
C     -1.3898      0.0572     -0.0114
H     -2.4836      0.1021     -0.0204
C     -0.6824      0.0925     -1.2088
H     -1.2194      0.1652     -2.1599
C      0.7075      0.0352     -1.1973
H      1.2641      0.0628     -2.1395
C      1.3899     -0.0572      0.0114
H      2.4836     -0.1022      0.0205
""")
close(file)

Dokoka.parseargs(false)

molecule1 = Dokoka.readmol("mock-build.xyz") 
molecule2 = Dokoka.readmol("mock-build.xyz")

Dokoka.goradial!(molecule1, molecule2, [0.0, 5.0, 0.0], [1.0, 0.0, 0.0], 0.0, 1)
rm("radial.xyz")

Dokoka.gorotational!(molecule1, molecule2, [0.0, 5.0, 0.0], [1.0, 0.0, 0.0], 0.0, 1)
rm("rotational.xyz")

Dokoka.goangular!(molecule1, molecule2, [0.0, 5.0, 0.0], [1.0, 0.0, 0.0], 0.0, 1)
rm("angular.xyz")

Dokoka.gorandom!(molecule1, molecule2, [0.0, 5.0, 0.0], [1.0, 0.0, 0.0], 0.0, 1)
rm("random.xyz")

rm("mock-build.xyz")
