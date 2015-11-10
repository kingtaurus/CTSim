function H = HKappa(E, EKappa, epsilonKappa)

H = heavisideLocal(E-EKappa, epsilonKappa);
