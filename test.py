# Git repo test
import math

def chi(z, rho):
    return math.log((math.sqrt(1 - 2 * rho * z + z * z) + z - rho) / (1 - rho))

def zfn(spotstrike_pow, nu, sigma0, log_moneyness):
    return (nu / sigma0) * spotstrike_pow * log_moneyness

def spotstrike_pow(f, K, beta):
    return pow(f * K, 0.5 * (1 - beta))

def black_iv(f, K, T, beta, nu, rho, sigma0):
    # Precompute
    log_moneyness = math.log(f / K)
    spotstrike_pow_ = spotstrike_pow(f, K, beta)
    z = zfn(spotstrike_pow_, nu, sigma0, log_moneyness)
    # Implied vol
    temp = 1 + pow(log_moneyness, 2) * pow(1 - beta, 2) / 24.0 + pow(log_moneyness, 4) * pow(1 - beta, 4) / 1920.0
    temp *= spotstrike_pow_
    out = (sigma0 / temp)
    out *= z / chi(z, rho)
    out *= 1 + ((pow(1 - beta, 2) / 24.0) * sigma0 * sigma0 / (spotstrike_pow_ * spotstrike_pow_)
                + rho * beta * nu * sigma0 / 4 / spotstrike_pow_
                + (2 - 3 * rho * rho) * nu * nu / 24.0) * T
    return out


def main():
    f = 100
    T = 1
    beta = 0.5
    nu = 0.1
    rho = -0.1
    sigma0 = 0.2

    IVs = []
    for K in range(60, 150):
        IVs.append(black_iv(f, K, T, beta, nu, rho, sigma0))

    for val in IVs:
        print(val)

if __name__ == "__main__":
    main()
