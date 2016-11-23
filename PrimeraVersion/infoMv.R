function (dp, x = NULL, y, w, penalty = NULL, norm2.tol = 1e-06) 
{
    type <- if (missing(y)) 
        "expected"
    else "observed"
    if (type == "observed") {
        if (!is.matrix(y)) 
            stop("y is not a matrix")
    }
    cp <- dp2cpMv(dp, "SN")
    d <- length(dp$alpha)
    d2 <- d * (d + 1)/2
    if (missing(w)) 
        w <- rep(1, max(NROW(cbind(x, y)), 1))
    if (any(w != round(w)) | any(w < 0)) 
        stop("weights must be non-negative integers")
    n <- length(w)
    nw <- sum(w)
    if (is.null(x)) {
        p <- 1
        xx <- sum.x <- nw
        x <- matrix(1, nrow = n, ncol = 1)
    }
    else {
        p <- NCOL(x)
        xx <- drop(t(x) %*% (w * x))
        sum.x <- drop(matrix(colSums(w * x)))
    }
    beta <- as.matrix(dp[[1]], p, d)
    Omega <- dp$Omega
    omega <- sqrt(diag(Omega))
    alpha <- dp$alpha
    eta <- alpha/omega
    Obar <- cov2cor(Omega)
    Obar.alpha <- as.vector(Obar %*% alpha)
    alpha.star <- sqrt(sum(alpha * Obar.alpha))
    if (alpha.star < 1e-04) {
        warning("information matrix of multivariate SN not computed near alpha=0")
        return(NULL)
    }
    c1 <- sqrt(2/pi)/sqrt(1 + alpha.star^2)
    c2 <- 1/(pi * sqrt(1 + 2 * alpha.star^2))
    D <- duplicationMatrix(d)
    i1 <- 1:prod(dim(beta))
    i2 <- max(i1) + 1:(d * (d + 1)/2)
    i3 <- max(i2) + 1:d
    O.inv <- pd.solve(Omega, silent = TRUE)
    if (type == "observed") {
        y0 <- y - x %*% beta
        S0 <- t(y0) %*% (w * y0)/nw
        y0.eta <- as.vector(y0 %*% eta)
        z1 <- zeta(1, y0.eta) * w
        z2 <- (-zeta(2, y0.eta) * w)
        S1 <- (O.inv %x% t(x)) %*% as.vector(w * y0) - (eta %x% 
            t(x)) %*% z1
        S2 <- (nw/2) * t(D) %*% ((O.inv %x% O.inv) %*% as.vector(S0 - 
            Omega))
        S3 <- t(y0) %*% z1
        score <- c(S1, S2, S3)
        u <- t(x) %*% z1
        U <- t(x) %*% (z2 * y0)
        V <- O.inv %*% (2 * S0 - Omega) %*% O.inv
        j11 <- O.inv %x% xx + outer(eta, eta) %x% (t(x) %*% (z2 * 
            x))
        j12 <- (O.inv %x% (t(x) %*% (w * y0) %*% O.inv)) %*% 
            D
        j13 <- diag(d) %x% u - eta %x% U
        j22 <- (nw/2) * t(D) %*% (O.inv %x% V) %*% D
        j23 <- matrix(0, d * (d + 1)/2, d)
        j33 <- t(y0) %*% (z2 * y0)
        uaA.coef <- NULL
    }
    else {
        Omega.eta <- omega * Obar.alpha
        mu.c <- Omega.eta/alpha.star^2
        Omega.c <- Omega - outer(Omega.eta, Omega.eta)/alpha.star^2
        alpha.bar <- alpha.star/sqrt(1 + 2 * alpha.star^2)
        ginvMills <- function(x, m = 0, s = 1) exp(-0.5 * ((x - 
            m)^2/s^2 - x^2) + log(zeta(1, x)) - log(s))
        fn.u <- function(x, sd, k) x^k * ginvMills(x, 0, sd)
        if (alpha.bar > 0) {
            err <- .Machine$double.eps^0.5
            u0 <- integrate(fn.u, -Inf, Inf, sd = alpha.bar, 
                k = 0, rel.tol = err)$value
            u1 <- integrate(fn.u, -Inf, Inf, sd = alpha.bar, 
                k = 1, rel.tol = err)$value
            u2 <- integrate(fn.u, -Inf, Inf, sd = alpha.bar, 
                k = 2, rel.tol = err)$value
        }
        else {
            u0 <- 2
            u1 <- u2 <- 0
        }
        a0 <- u0
        a1 <- u1 * mu.c
        A2 <- u2 * outer(mu.c, mu.c) + u0 * Omega.c
        A1 <- (c1 * (diag(d) - outer(eta, eta) %*% Omega/(1 + 
            alpha.star^2)) - c2 * outer(eta, a1))
        j11 <- (O.inv + c2 * a0 * outer(eta, eta)) %x% xx
        j12 <- c1 * (O.inv %x% outer(sum.x, eta)) %*% D
        j13 <- A1 %x% sum.x
        j22 <- 0.5 * nw * t(D) %*% (O.inv %x% O.inv) %*% D
        j23 <- matrix(0, d * (d + 1)/2, d)
        j33 <- nw * c2 * A2
        uaA.coef <- list(u0 = u0, u1 = u1, u2 = u2, a1 = a1, 
            A1 = A1, A2 = A2)
        score <- NULL
    }
    I.theta <- rbind(cbind(j11, j12, j13), cbind(t(j12), j22, 
        j23), cbind(t(j13), t(j23), j33))
    if (!is.null(penalty)) {
        penalty.fn <- if (is.null(penalty)) 
            NULL
        else get(penalty, inherits = TRUE)
        penalty.theta <- function(theta23, penalty, d) {
            vOmega <- theta23[1:(d * (d + 1)/2)]
            eta <- theta23[(d * (d + 1)/2) + (1:d)]
            Omega <- vech2mat(vOmega)
            alpha <- eta * sqrt(diag(Omega))
            penalty(list(alpha = alpha, Omega = Omega))
        }
        i23 <- c(i2, i3)
        theta23 <- c(Omega[lower.tri(Omega, TRUE)], eta)
        score[i23] <- (score[i23] - numDeriv::grad(penalty.theta, 
            theta23, penalty = penalty.fn, d = d))
        jQ <- numDeriv::hessian(penalty.theta, theta23, penalty = penalty.fn, 
            d = d)
        I.theta[i23, i23] <- I.theta[i23, i23] + jQ
    }
    I.theta <- force.symmetry(I.theta, tol = 1000)
    inv_I.theta <- pd.solve(I.theta, silent = TRUE)
    if (is.null(inv_I.theta)) {
        warning("numerically unstable information matrix")
        return(NULL)
    }
    if (type == "observed") {
        score.norm2 <- sum(score * as.vector(inv_I.theta %*% 
            score))
        if (score.norm2/d > norm2.tol) 
            warning("'dp' does not seem to be at MLE")
    }
    D32 <- matrix(0, d, d2)
    tmp32 <- matrix(0, d^2, d^2)
    for (i in 1:d) {
        Eii <- matrix(0, d, d)
        Eii[i, i] <- 1
        tmp32 <- tmp32 + Eii %x% Eii
    }
    D32 <- (-0.5) * (t(eta) %x% diag(1/omega^2, d, d)) %*% tmp32 %*% 
        D
    Dlow <- cbind(matrix(0, d, d * p), D32, diag(1/omega, d, 
        d))
    Dtheta.dp <- rbind(cbind(diag(d * p + d2), matrix(0, d * 
        p + d2, d)), Dlow)
    I.dp <- t(Dtheta.dp) %*% I.theta %*% Dtheta.dp
    I.dp <- force.symmetry(I.dp, tol = 1000)
    Sigma <- cp$var.cov
    sigma <- sqrt(diag(Sigma))
    Sigma.inv <- pd.solve(Sigma)
    mu0 <- c1 * omega * Obar.alpha
    beta0.sq <- as.vector(t(mu0) %*% Sigma.inv %*% mu0)
    beta0 <- sqrt(beta0.sq)
    q1 <- 1/(c1 * (1 + beta0.sq))
    q2 <- 0.5 * q1 * (2 * c1 - q1)
    Dplus <- pd.solve(t(D) %*% D) %*% t(D)
    D23 <- Dplus %*% (diag(d) %x% mu0 + mu0 %x% diag(d))
    a <- as.vector(Sigma.inv %*% mu0)
    D32 <- t(-a) %x% (q1 * Sigma.inv - q1 * q2 * outer(a, a)) %*% 
        D
    D33 <- q1 * Sigma.inv - 2 * q1 * q2 * outer(a, a)
    one00 <- c(1, rep(0, p - 1))
    Dtheta.psi <- rbind(cbind(diag(p * d), matrix(0, p * d, d2), 
        -diag(d) %x% one00), cbind(matrix(0, d2, p * d), diag(d2), 
        D23), cbind(matrix(0, d, p * d), D32, D33))
    mu0. <- mu0/(sigma * beta0)
    D32. <- matrix(0, d, d2)
    for (i in 1:d) {
        Eii <- matrix(0, d, d)
        Eii[i, i] <- 1
        D32. <- D32. + (1/sigma[i]) * ((t(mu0.) %*% Eii) %x% 
            Eii) %*% D
    }
    D32. <- 0.5 * beta0 * D32.
    D33. <- (2/(4 - pi)) * diag(sigma/mu0.^2, d, d)/(3 * beta0.sq)
    Dpsi.cp <- rbind(cbind(diag(p * d + d2), matrix(0, p * d + 
        d2, d)), cbind(matrix(0, d, p * d), D32., D33.))
    jacob <- Dtheta.psi %*% Dpsi.cp
    I.cp <- t(jacob) %*% I.theta %*% jacob
    I.cp <- if (any(is.na(I.cp))) 
        NULL
    else force.symmetry(I.cp)
    asyvar.dp <- pd.solve(I.dp, silent = TRUE)
    if (is.null(asyvar.dp)) 
        se.dp <- list(NULL)
    else {
        diags.dp <- sqrt(diag(asyvar.dp))
        se.beta <- matrix(diags.dp[1:(p * d)], p, d)
        se.diagOmega <- diags.dp[p * d + d2 + 1 - rev(cumsum(1:d))]
        se.alpha <- diags.dp[p * d + d2 + (1:d)]
        se.dp <- list(beta = se.beta, diagOmega = se.diagOmega, 
            alpha = se.alpha)
    }
    asyvar.cp <- pd.solve(I.cp, silent = TRUE)
    if (is.null(asyvar.cp)) 
        se.cp <- list(NULL)
    else {
        diags.cp <- sqrt(diag(asyvar.cp))
        se.beta <- matrix(diags.cp[1:(p * d)], p, d)
        se.diagSigma <- diags.cp[p * d + d2 + 1 - rev(cumsum(1:d))]
        se.gamma1 <- diags.cp[p * d + d2 + (1:d)]
        se.cp <- list(beta = se.beta, var = se.diagSigma, gamma1 = se.gamma1)
    }
    aux <- list(info.theta = I.theta, score.theta = score, Dtheta.dp = Dtheta.dp, 
        Dpsi.cp = Dpsi.cp, Dtheta.psi = Dtheta.psi, uaA.coef = uaA.coef)
    list(dp = dp, cp = cp, type = type, info.dp = I.dp, info.cp = I.cp, 
        asyvar.dp = asyvar.dp, asyvar.cp = asyvar.cp, se.dp = se.dp, 
        se.cp = se.cp, aux = aux)
}
