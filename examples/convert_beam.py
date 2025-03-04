import math
import numpy as np

def geodetic_to_ecef(lat_deg, lon_deg, h):
    """
    Convert geodetic coordinates (lat, lon in degrees, elevation in feet)
    to ECEF (Earth-Centered, Earth-Fixed) Cartesian coordinates in meters.
    """
    # WGS84 ellipsoid constants
    a = 6378137.0  # semi-major axis in meters
    f = 1 / 298.257223563  # flattening
    e2 = 2 * f - f**2      # eccentricity squared

    # Convert latitude and longitude to radians
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)

    # Prime vertical radius of curvature
    sin_lat = math.sin(lat)
    N = a / math.sqrt(1 - e2 * sin_lat**2)

    # Compute ECEF coordinates
    x = (N + h) * math.cos(lat) * math.cos(lon)
    y = (N + h) * math.cos(lat) * math.sin(lon)
    z = ((1 - e2) * N + h) * sin_lat

    return np.array([x, y, z])

def kabsch(P, Q):
    """
    Computes the optimal rotation matrix R and translation vector T
    that maps points P (Nx3) onto points Q (Nx3) using the Kabsch algorithm.
    The mapping is: Q = R * P + T
    """
    # Compute centroids
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)

    # Center the points
    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q

    # Compute covariance matrix
    H = P_centered.T @ Q_centered

    # Singular Value Decomposition
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T

    # Ensure a proper rotation (no reflection)
    if np.linalg.det(R) < 0:
        Vt[2, :] *= -1
        R = Vt.T @ U.T

    # Compute translation
    T = centroid_Q - R @ centroid_P
    return R, T

if __name__ == "__main__":
    # Reference input: three GPS points (lat, lon, elevation in feet)
    gps_points = [
        (41.83743939, -88.2694776399, 740.6*0.3048), # microboone
        (47.82027, -92.24141, -248.3992), # minos FD
        (48.37859633, -92.83128217, 348.638) # nova FD
    ]

    # Convert GPS coordinates to ECEF coordinates
    ecef_points = np.array([geodetic_to_ecef(lat, lon, elev) for lat, lon, elev in gps_points])

    # Reference corresponding points in the new (target) Cartesian system
    # (These should be provided or measured in your new coordinate system.)
    target_points = np.array([
        [53., 76., 679.],
        [0, 0, 735340.],
        [11037.296, -4162.557, 810422.32]
    ])

    # Compute the mapping using the Kabsch algorithm
    R, T = kabsch(ecef_points, target_points)

    print("Rotation Matrix R:")
    print(R)
    print("\nTranslation Vector T:")
    print(T)

    # Now, for a new GPS coordinate, compute its ECEF coordinates and map to the target system.
    lp1 = (47.016657, -91.647625, 179) # lake point 1
    lp2 = (47.004683, -91.665225, 179) # lake point 2
    th  = (47.0214787, -91.664325, 203) # two harbors
    lp3 = (46.980478, -91.610088, 179) # lake point 3 - best
    new_gps = lp3

    new_ecef = geodetic_to_ecef(*new_gps)
    print("\nNew GPS Coordinate (ECEF system):")
    print(new_ecef)

    new_target = R @ new_ecef + T
    print("\nNew Coordinate (ECEF mapped to target system) (cm):")
    print(100.*new_target)
