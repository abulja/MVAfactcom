#load data
com = read.table(file="BusinessCompetitiveness.txt")

#define x
x=as.matrix(com)

#define number of row and col
n=nrow(x)
p=ncol(x)

#define variable names
colnames(x) = c("QES","QMSE","FLTA","FDI","NCA","PPS","CI","QSRI","CSR&D","UICR&D","GPATP","ASE","PCTP")

#correlation matrix
r = cor(x)
m = r

for (i in 1:ncol(r)) {
    m[i, i] = r[i, i] - 1
}

psi = matrix(1,13,13)
for (i in 1:13) {
    psi[i, i] = 1 - max(abs(m[, i]))
}

#spectral decomposition
eig = eigen(r - psi)
ee  = eig$values[1:2]
vv  = eig$vectors[, 1:2]
vv  = t(t(vv[, 1:2]) * sign(vv[2, 1:2]))
q1  = sqrt(ee[1]) * vv[, 1]
q2  = sqrt(ee[2]) * vv[, 2]
q   = cbind(q1, q2)

#plot
plot(q, type = "n", xlab = "First Factor", ylab = "Second Factor", main = "Business competitiveness in small open economies", 
    xlim = c(-1.8, 1), ylim = c(-0.7, 0.7), cex.lab = 1.4, cex.axis = 1.4, cex.main = 1.8)
text(q, colnames(x), cex = 1.2, xpd = NA)
abline(v = 0)
abline(h = 0)

#Maximum Likelihood Factor Analysis without varimax rotation factanal performs
mlm  = factanal(xs, 3, rotation = "none", covmat = r)
load = mlm$loadings                           #estimated factor loadings
ld   = cbind(load[, 1], load[, 2], load[, 3]) #the estimated factor loadings matrix
com  = diag(ld %*% t(ld))                     #communalities are calculated
psi  = diag(r) - diag(ld %*% t(ld))         #specific variances are calculated
tbl  = cbind(load[, 1], load[, 2], load[, 3], com, psi)

dev.new()
par(mfcol = c(2, 2))

#plot first factor against second
plot(load[, 1], load[, 2], type = "n", xlab = "x", ylab = "y", main = "Factors21 - theta21", 
    font.main = 1, cex.lab = 1.1, cex.axis = 1.1, cex.main = 1.4, ylim = c(-0.6, 0.6))
text(load[, 1], load[, 2], colnames(x), cex = 1.1)
abline(h = 0, v = 0)

#plot first factor against third
plot(load[, 1], load[, 3], type = "n", xlab = "x", ylab = "y", main = "Factors31 - theta31", 
    font.main = 1, cex.lab = 1.1, cex.axis = 1.1, cex.main = 1.4, ylim = c(-0.4, 0.4))
text(load[, 1], load[, 3], colnames(x), cex = 1.1)
abline(h = 0, v = 0)

#plot second factor against third
plot(load[, 2], load[, 3], type = "n", xlab = "x", ylab = "y", main = "Factors32 - theta32", 
    font.main = 1, cex.lab = 1.1, cex.axis = 1.1, cex.main = 1.4, xlim = c(-0.6, 0.6), 
    ylim = c(-0.4, 0.4))
text(load[, 2], load[, 3], colnames(x), cex = 1.1)
abline(h = 0, v = 0)

#Maximum Likelihood Factor Analysis after varimax rotation
var  = varimax(ld)                            #rotates the factor loadings matrix
load = var$loadings                           #estimated factor loadings after varimax
vl   = cbind(load[, 1], load[, 2], load[, 3])
com  = diag(vl %*% t(vl))                     #communalities are calculated
psi  = diag(r) - diag(vl %*% t(vl))         #specific variances are calculated
tbl  = cbind(load[, 1], load[, 2], load[, 3], com, psi)

dev.new()
par(mfcol = c(2, 2))

#plot first factor against second
plot(load[, 1], load[, 2], type = "n", xlab = "x", ylab = "y", main = "Factors21 - theta21", 
    font.main = 1, cex.lab = 1.1, cex.axis = 1.1, cex.main = 1.4, xlim = c(-1, 1))
text(load[, 1], load[, 2], colnames(x), cex = 1.1)
abline(h = 0, v = 0)

#plot first factor against third
plot(load[, 1], load[, 3], type = "n", xlab = "x", ylab = "y", main = "Factors31 - theta31", 
    font.main = 1, cex.lab = 1.1, cex.axis = 1.1, cex.main = 1.4, xlim = c(-1, 1))
text(load[, 1], load[, 3], colnames(x), cex = 1.1)
abline(h = 0, v = 0)

#plot second factor against third
plot(load[, 2], load[, 3], type = "n", xlab = "x", ylab = "y", main = "Factors32 - theta32", 
    font.main = 1, cex.lab = 1.1, cex.axis = 1.1, cex.main = 1.4, xlim = c(-1, 1))
text(load[, 2], load[, 3], colnames(x), cex = 1.1)
abline(h = 0, v = 0)

#Principal Component Method after varimax rotation spectral decomposition
e      = eigen(r)
eigval = e$values[1:3]
eigvec = e$vectors[, 1:3]
E      = matrix(eigval, nrow(r), ncol = 3, byrow = T)
Q      = sqrt(E) * eigvec                     #the estimated factor loadings matrix
pcm    = varimax(Q)                           #rotates the factor loadings matrix
load   = pcm$loadings                         #estimated factor loadings after varimax
ld     = cbind(load[, 1], load[, 2], load[, 3])
com    = diag(ld %*% t(ld))                   #communalities are calculated
psi    = diag(r) - diag(ld %*% t(ld))       #specific variances are calculated
tbl    = cbind(load[, 1], load[, 2], load[, 3], com, psi)

dev.new()
par(mfcol = c(2, 2))

#plot first factor against second
plot(load[, 1], load[, 2], type = "n", xlab = "x", ylab = "y", main = "Factors21 - theta21", 
    font.main = 1, cex.lab = 1.1, cex.axis = 1.1, cex.main = 1.4)
text(load[, 1], load[, 2], colnames(x), cex = 1.1)
abline(h = 0, v = 0)

#plot first factor against third
plot(load[, 1], load[, 3], type = "n", xlab = "x", ylab = "y", main = "Factors31 - theta31", 
    font.main = 1, cex.lab = 1.1, cex.axis = 1.1, cex.main = 1.4)
text(load[, 1], load[, 3], colnames(x), cex = 1.1)
abline(h = 0, v = 0)

#plot second factor against third
plot(load[, 2], load[, 3], type = "n", xlab = "x", ylab = "y", main = "Factors32 - theta32", 
    font.main = 1, cex.lab = 1.1, cex.axis = 1.1, cex.main = 1.4)
text(load[, 2], load[, 3], colnames(x), cex = 1.1)
abline(h = 0, v = 0)

#Principal Factor Method after varimax rotation inverse of the correlation matrix
f      = solve(r)
psiini = diag(1/f[row(f) == col(f)]) #preliminary estimate of psi
psi    = psiini
for (i in 1:10) {
    ee     = eigen(r - psi)
    eigval = ee$values[1:3]
    eigvec = ee$vectors[, 1:3]
    EE     = matrix(eigval, nrow(r), ncol = 3, byrow = T)
    QQ     = sqrt(EE) * eigvec
    psiold = psi
    psi    = diag(as.vector(1 - t(colSums(t(QQ * QQ)))))
    i      = i + 1
    z      = psi - psiold
    convergence = z[row(z) == col(z)]
}
pfm  = varimax(QQ)                             #rotates the factor loadings matrix
load = pfm$loadings                            #estimated factor loadings after varimax
ld   = cbind(load[, 1], load[, 2], load[, 3])
com  = diag(ld %*% t(ld))                      #communalities are calculated
psi  = diag(r) - diag(ld %*% t(ld))          #specific variances are calculated
tbl  = cbind(load[, 1], load[, 2], load[, 3], com, psi)

dev.new()
par(mfcol = c(2, 2))

#plot first factor against second
plot(load[, 1], load[, 2], type = "n", xlab = "x", ylab = "y", main = "Factors21 - theta21", 
    font.main = 1, cex.lab = 1.1, cex.axis = 1.1, cex.main = 1.4)
text(load[, 1], load[, 2], colnames(x), cex = 1.1)
abline(h = 0, v = 0)

#plot first factor against third
plot(load[, 1], load[, 3], type = "n", xlab = "x", ylab = "y", main = "Factors31 - theta31", 
    font.main = 1, cex.lab = 1.1, cex.axis = 1.1, cex.main = 1.4, ylim = c(-1, 1))
text(load[, 1], load[, 3], colnames(x), cex = 1.1)
abline(h = 0, v = 0)

#plot second factor against third
plot(load[, 2], load[, 3], type = "n", xlab = "x", ylab = "y", main = "Factors32 - theta32", 
    font.main = 1, cex.lab = 1.1, cex.axis = 1.1, cex.main = 1.4, ylim = c(-1, 1))
text(load[, 2], load[, 3], colnames(x), cex = 1.1)
abline(h = 0, v = 0)

    


