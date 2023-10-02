import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
import math

THRUST = 8000  # N
PRESSURE = 50 * 100000  # Pa

L_STAR = 0.45

CONTRACTION_RATIO = 6
EXPANSION_RATIO = 10

SEALEVEL_PRESSURE = 1 * 100000


COMBUSION_TEMPERATURE = 3540  # K
Y = 1.220  # Somwhere between 1.2 and 1.4 i think?

# Constants
PI = 3.1415
G = 9.81
R = 8.31446261815324


THEORETICAL_CHARACHTERISTIC_VELOCITY = (
    (G * Y * R * COMBUSION_TEMPERATURE) ** (1 / 2)
) / (Y * (((2) / (Y + 1)) ** ((Y + 1) / (Y - 1))) ** (1 / 2))
DESIGN_CHARACHTERISTIC_VELOCITY = THEORETICAL_CHARACHTERISTIC_VELOCITY * 0.975


def thrust_coefficient(ambient_pressure, exit_pressure):
    return (
        ((2 * Y**2) / (Y - 1))
        * (((2) / (Y + 1)) ** ((Y + 1) / (Y - 1)))
        * ((1 - ((exit_pressure) / (PRESSURE))) ** ((Y - 1) / (Y)))
    ) ** (1 / 2) + EXPANSION_RATIO * ((exit_pressure - ambient_pressure) / (PRESSURE))


THEORETICAL_SEALEVEL_THRUST_COEFFICIENT = thrust_coefficient(
    SEALEVEL_PRESSURE, PRESSURE * (1 / EXPANSION_RATIO)
)
DESIGN_SEALEVEL_THRUST_COEFFICIENT = THEORETICAL_SEALEVEL_THRUST_COEFFICIENT * 0.98

DESIGN_SEALEVEL_ISP = (
    DESIGN_CHARACHTERISTIC_VELOCITY * DESIGN_SEALEVEL_THRUST_COEFFICIENT
) / (G)

print()
print(
    f"Sealevel:\nThrust Coefficient: {DESIGN_SEALEVEL_THRUST_COEFFICIENT:.3f}\nc: {DESIGN_CHARACHTERISTIC_VELOCITY:.2f}m/s\nIsp: {DESIGN_SEALEVEL_ISP}s"
)
print()

THROAT_AREA = THRUST / (DESIGN_SEALEVEL_THRUST_COEFFICIENT * PRESSURE)
THROAT_RADIUS = ((THROAT_AREA) / (PI)) ** (1 / 2)

EXIT_RADIUS = THROAT_RADIUS * (EXPANSION_RATIO) ** (1 / 2)

CHAMBER_VOLUME = THROAT_AREA * L_STAR

print(
    f"Inner dimensions:\nThroat area: {THROAT_AREA:.6f}m^2\nThroat radius: {THROAT_RADIUS:.6f}m\nExit radius: {EXIT_RADIUS:.6f}m\nChamber volume: {CHAMBER_VOLUME:.6f}m^3"
)
print()

CONVERGING_THROAT_SECTION_RADIUS = 1.5 * THROAT_RADIUS
DIVERGING_THROAT_SECTION_RADIUS = 0.382 * THROAT_RADIUS

CHAMBER_RADIUS = THROAT_RADIUS * (CONTRACTION_RATIO) ** (1 / 2)

CONVERGENT_CONE_LENGTH = (
    THROAT_RADIUS * ((CONTRACTION_RATIO) ** (1 / 2) - 1)
    + CONVERGING_THROAT_SECTION_RADIUS * ((1 / math.cos(math.radians(20))) - 1)
) / (math.tan(math.radians(20)))

CONVERGENT_CONE_VOLUME = (
    (PI / 3)
    * CONVERGENT_CONE_LENGTH
    * (CHAMBER_RADIUS**2 + THROAT_RADIUS**2 + CHAMBER_RADIUS * THROAT_RADIUS)
)

CYLINDRICAL_CHAMBER_LENGTH = (CHAMBER_VOLUME - CONVERGENT_CONE_VOLUME) / (
    CONTRACTION_RATIO * THROAT_AREA
)

NOZZLE_LENGTH = (
    0.8
    * (
        THROAT_RADIUS * ((EXPANSION_RATIO) ** (1 / 2) - 1)
        + DIVERGING_THROAT_SECTION_RADIUS * (1 / math.cos(math.radians(15)) - 1)
    )
    / (math.tan(math.radians(15)))
)

plt.plot(
    [
        -(CYLINDRICAL_CHAMBER_LENGTH + CONVERGENT_CONE_LENGTH),
        -(CYLINDRICAL_CHAMBER_LENGTH + CONVERGENT_CONE_LENGTH),
        -CONVERGENT_CONE_LENGTH,
    ],
    [0, CHAMBER_RADIUS, CHAMBER_RADIUS],
)

plt.plot([-CONVERGENT_CONE_LENGTH, 0], [CHAMBER_RADIUS, THROAT_RADIUS], "-.")

plt.plot(
    -(CYLINDRICAL_CHAMBER_LENGTH + CONVERGENT_CONE_LENGTH),
    0,
    "ro",
)

plt.plot(-(CYLINDRICAL_CHAMBER_LENGTH + CONVERGENT_CONE_LENGTH), CHAMBER_RADIUS, "ro")
plt.plot(-CONVERGENT_CONE_LENGTH, CHAMBER_RADIUS, "ro")

plt.plot(0, THROAT_RADIUS, "ro")

plt.plot(
    NOZZLE_LENGTH,
    EXIT_RADIUS,
    "ro",
)

N_ANGLE = 27.4
E_ANGLE = 9.8

DIVERGING_THROAT_LENGTH = DIVERGING_THROAT_SECTION_RADIUS * math.sin(
    math.radians(N_ANGLE)
)

DIVERGING_THROAT_END_RADIUS = THROAT_RADIUS + DIVERGING_THROAT_SECTION_RADIUS * (
    1 - math.cos(math.radians(N_ANGLE))
)


def equations(p):
    a, b, c, d = p
    return (
        a * DIVERGING_THROAT_LENGTH**3
        + b * DIVERGING_THROAT_LENGTH**2
        + c * DIVERGING_THROAT_LENGTH
        + d
        - DIVERGING_THROAT_END_RADIUS,
        a * NOZZLE_LENGTH**3
        + b * NOZZLE_LENGTH**2
        + c * NOZZLE_LENGTH
        + d
        - EXIT_RADIUS,
        3 * a * DIVERGING_THROAT_LENGTH**2
        + 2 * b * DIVERGING_THROAT_LENGTH**1
        + c
        - math.tan(math.radians(N_ANGLE)),
        3 * a * NOZZLE_LENGTH**2
        + 2 * b * NOZZLE_LENGTH**1
        + c
        - math.tan(math.radians(E_ANGLE)),
    )


a, b, c, d = fsolve(equations, (1, 1, 1, 1))


def Mach_Area_Relation(A_over_At, sup):
    def equations(p):
        M = p
        return (1 / M) * (
            ((2 / (Y + 1)) * (1 + ((Y - 1) / 2) * M**2)) ** ((Y + 1) / (2 * (Y - 1)))
        ) - A_over_At

    M = fsolve(equations, (0.1 + 100 * int(sup)))
    print(M)
    return M


def Mach(x):
    return Mach_Area_Relation((PI * (radius(x) ** 2)) / THROAT_AREA, x > 0)


def Pressure(x):
    return PRESSURE * ((1 + (((Y - 1) / (2)) * (Mach(x) ** 2))) ** ((-Y) / (Y - 1)))


def Temp(x):
    return COMBUSION_TEMPERATURE * ((1 + (((Y - 1) / (2)) * (Mach(x) ** 2))) ** (-1))


def Density(x):
    return 5 * ((1 + (((Y - 1) / (2)) * (Mach(x) ** 2))) ** ((-Y) / (Y - 1)))


plt.plot(
    DIVERGING_THROAT_SECTION_RADIUS * math.sin(math.radians(N_ANGLE)),
    DIVERGING_THROAT_END_RADIUS,
    "ro",
)


def radius(x):
    if x <= 0 and not (
        (CONVERGING_THROAT_SECTION_RADIUS**2 > x**2)
        and (
            (
                CHAMBER_RADIUS
                >= THROAT_RADIUS
                + CONVERGING_THROAT_SECTION_RADIUS
                - (((CONVERGING_THROAT_SECTION_RADIUS**2) - (x**2)) ** (1 / 2))
            )
        )
    ):
        return CHAMBER_RADIUS
    elif x <= 0:
        return (THROAT_RADIUS + CONVERGING_THROAT_SECTION_RADIUS) - (
            (CONVERGING_THROAT_SECTION_RADIUS**2) - (x**2)
        ) ** (1 / 2)
    elif x <= DIVERGING_THROAT_LENGTH:
        return (THROAT_RADIUS + DIVERGING_THROAT_SECTION_RADIUS) - (
            (DIVERGING_THROAT_SECTION_RADIUS**2) - (x**2)
        ) ** (1 / 2)
    else:
        return a * x**3 + b * x**2 + c * x + d


x = np.linspace(
    -(CYLINDRICAL_CHAMBER_LENGTH + CONVERGENT_CONE_LENGTH), NOZZLE_LENGTH, 10000
)
y = np.array([radius(a) for a in x])
plt.plot(x, y, color="black", label="Contour")

mach = np.array([Mach(a) for a in x])
plt.plot(x, mach * 0.01, color="orange", label="Mach")

pressure = np.array([Pressure(a) for a in x])
plt.plot(x, pressure * (1 / 100000) * 0.0012, color="blue", label="Pressure")

temp = np.array([Temp(a) for a in x])
plt.plot(x, temp * (1 / 1000) * 0.01, color="red", label="Temp")

density = np.array([Density(a) for a in x])
plt.plot(x, density * 0.01, color="green", label="Density")

Hg = np.array([(Mach(a) * Density(a)) ** (4 / 5) for a in x])
plt.plot(x, Hg * 0.01, color="purple", label="Hg")

plt.legend(loc="upper left")

plt.xlim(
    -(
        CYLINDRICAL_CHAMBER_LENGTH
        + CONVERGENT_CONE_LENGTH
        + CONVERGENT_CONE_LENGTH * 0.1
    ),
    NOZZLE_LENGTH + CONVERGENT_CONE_LENGTH * 0.1,
)
plt.ylim(0)
plt.gca().set_aspect("equal", adjustable="box")

plt.show()
