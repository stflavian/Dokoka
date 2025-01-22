using Dokoka

Dokoka.parseargs(false)
#Dokoka.goradial!(Dokoka.readmol("../benzene.xyz"), Dokoka.readmol("../benzene.xyz"), [0.0, 5.0, 0.0], [1.0, 0.0, 0.0], 0.0, 1)
#Dokoka.gorotational!(molecule1, molecule2, [0.0, 5.0, 0.0], [1.0, 0.0, 0.0], 0.0, 1)
#Dokoka.goangular!(molecule1, molecule2, [0.0, 5.0, 0.0], [1.0, 0.0, 0.0], 0.0, 1)
#Dokoka.gorandom!(molecule1, molecule2, [0.0, 5.0, 0.0], [1.0, 0.0, 0.0], 0.0, 1)

