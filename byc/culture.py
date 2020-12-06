import numpy as np

def starting_OD(target_OD, duration_hrs, doubling_time_hrs=1.5):
    """
    Return start_OD (A_naught) that cells should be diluted to so that they
    will arrive at target_OD (A(t)) after duration_hrs given a 
    doubling time of doubling_time_hrs (defaults to 1.5)

    By A(t) = A_naught*e^kt,

        k = ln(2)/t 
        where t = doubling time (t when 2*A_naught = A_naught*e^kt)

    And,

        A_naught = A(t)/e^kt
    """
    k = np.log(2)/doubling_time_hrs
    start_OD = target_OD/np.exp(k*duration_hrs)

    return start_OD