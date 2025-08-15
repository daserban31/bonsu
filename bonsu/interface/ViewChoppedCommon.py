"""
Created on Fri Aug 15 12:17:24 2025.

@author: davids_work
"""
################################################################################
#                                 Imports                                      #
################################################################################
import numpy


################################################################################
#                                Constants                                     #
################################################################################


################################################################################
#                                 Methods                                      #
################################################################################
def normalise_vector(input_vector: numpy.ndarray) -> numpy.ndarray:
    """
    Return the normalised <input_vector>.

    Parameters
    ----------
    input_vector: numpy.ndarray
        The vector to be normalised.

    Return
    ------
    numpy.ndarray
        The vector normalised.
    """
    return numpy.array(input_vector)/numpy.linalg.norm(input_vector)


def skew_vector(vector: numpy.ndarray) -> numpy.ndarray:
    """
    Return the skewed (anti-symmetric) matrix from <vector>.

    Parameters
    ----------
    vector: numpy.ndarray
        Vector from which to calculate its anti-symmetric matrix, necessary
        for the calculation of the rotation matrix around it.

    Return
    ------
    numpy.ndarray
        The anti-symmetric matrix
    """
    return numpy.array((( 0,         -vector[2],  vector[1]), # noqa
                        ( vector[2],          0, -vector[0]), # noqa
                        (-vector[1],  vector[0],          0)))


def calc_rot_matrix(angle_rad: float, axis: numpy.ndarray) -> numpy.ndarray:
    """
    Return the rotation matrix of <angle_rad> around <axis>.

    Parameters
    ----------
    angle_rad: float
        Angle by which to rotate (in rad).
    axis: numpy.ndarray
        Vector around which the rotation is done.

    Return
    ------
    numpy.ndarray
        The rotation matrix
    """
    axis = normalise_vector(axis)
    return (numpy.cos(angle_rad)*numpy.identity(3) +
            numpy.sin(angle_rad)*skew_vector(axis) +
            (1 - numpy.cos(angle_rad))*numpy.outer(axis, axis))


def rotate_vector(vector: numpy.ndarray, normal: numpy.ndarray,
                  angle_deg: float) -> numpy.ndarray:
    """
    Return the <vector> rotated by <angle_deg> around <normal>.

    Parameters
    ----------
    vector: numpy.ndarray
        Initial vector to be rotated.
    normal: numpy.ndarray
        Axis around which to rotate, normal to the plane of rotation.
    angle_deg: float
        Angle by which to rotate (in degrees).

    Return
    ------
    numpy.ndarray
        The rotated vector.
    """
    return numpy.dot(calc_rot_matrix(numpy.deg2rad(angle_deg), normal),
                     vector)


def rotate_zhat_by_theta_psi(theta_deg: float,
                             psi_deg: float) -> numpy.ndarray:
    """Rotate z-hat by polar and azimuthal angles, respectively."""
    return rotate_vector(rotate_vector((0, 0, 1), (0, 1, 0), theta_deg),
                         (0, 0, 1), psi_deg)


################################################################################
#                                 Classes                                      #
################################################################################


################################################################################
#                         Dunder Name Dunder Main                              #
################################################################################
if __name__ == '__main__':
    pass
