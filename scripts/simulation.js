/* eslint-disable no-unused-vars */
/* global SceneUtils */
// setupGUI - turn it into a module and delete the above line
/** Simulation of a Mars lander.
 *
 *
 *
 * 6DOF Entry Trajectory Code with Geometry & Aerodynamics
 * @author Patryk Szczepanski
 */

// ///////////////////////////////       /////////////////////////////////////////////
// /////////////////////////////// MODEL /////////////////////////////////////////////
// ///////////////////////////////       /////////////////////////////////////////////

// /////////////////////
// //LANDER (MODULE 1)//
// /////////////////////

// define lander object with defineMesh method using lander parameters
const lander = {
  mass: 980, // * entry mass of the lander (kg)
  diameter: 2.65, // * diameter of the lander [m]
  get rBody() { return this.diameter / 2; }, // * body radius [m]
  rNose: 0.88, // * nose radius [m]
  get rShoulder() { return this.rNose / 10; }, // * shoulder radius [m]
  sphereConeAngleDeg: 70, // * angle between conic section and body x-axis [m]
  get theta() { return this.sphereConeAngleDeg * Math.PI / 180; }, // * same angle in radians
  get refArea() { return Math.PI * this.rBody * this.rBody; }, // reference surface area [m^2]
  COMx: 0.68, // payload x coordinate of COM
  get COMpayload() { return [-this.COMx, 0, 0]; }, // payload (cylinder) centre of mass [m]
  get rPayload() { return this.rBody / 1.6; }, // payload container radius [m]
  heightPayload: 1.4, // payload container height [m]
  // define heatshield properties
  hsThickness: 0.07, // [m]
  hsDensity: 265, // [kg/m3]
  // mesh of the lander
  mesh: {
    scale: 100000, // variable defining scale 1:value between real results and graphics dimensions
    landerMagnification: 3000,
    nop: 20, // number of mesh points in two directions of the mesh surface
    x: null, // x, y, z matrices defining mesh coordinates used for physics calculation
    y: null,
    z: null,
    xLath: null, // xLath and rLath defining the curve used in graphics module for the mesh
    rLath: null,
  },
  /**
   * Returns 2-D grid coordinates based on the coordinates contained in vectors x and y. X is a matrix where
   * each row is a copy of x, and Y is a matrix where each column is a copy of y. (same as meshgrid in matlab)
   *
   * @param {array} x - x-coordinates of points, specified as a vector
   * @param {array} y - y-coordinates of points, specified as a vector
   *
   * @return {array} - the array with matrices X and Y
   */
  meshgrid(x, y) {
    // initializes final matrices
    const X = [];
    const Y = [];

    // defines lengths of the vectors
    const xLength = x.length;
    const yLength = y.length;

    // creates matrices X and Y
    for (let i = 0; i < yLength; i++) {
      X.push(x); // adds rows of X
      const temp = []; // creates an temporary array to store rows of Y matrix

      for (let j = 0; j < xLength; j++) {
        temp.push(y[i]); // generates rows of Y
      }
      Y.push(temp); // adds rows of Y
    }
    return [X, Y];
  },
  /**
   * Calculates mesh nodes of the vehicle in the body frame x, y, z coordinates.
   * Also defines curve used in graphics module to create the mesh.
   *
   * @param {object} this - lander object
   *
   * @return {array} - array with x, y, z coordinates of the mesh and xLath, rLath curve coordinates.
   */
  defineMesh() {
    // defines variables from the lander object
    const { rBody, rNose, rShoulder, theta } = this;

    // define dimensions for nose, frustum and shoulder
    // nose
    // const noseArea = Math.PI * (rNose * Math.cos(theta)) * (rNose * Math.cos(theta));// area of projected nose [m^2] //UNUSED
    const noseH = rNose * (1 - Math.sin(theta)); // X-coordinate where nose edge meets frustum [m]
    const noseL = rNose * Math.cos(theta); // R-coordinate where nose edge meets frustum [m]
    // const noseCap = 2 * Math.PI * rNose * noseH; // area of a spherical cap [m^2] //UNUSED

    // shoulder height and length
    const shoulderH = rShoulder * (1 - Math.cos(theta));
    const shoulderL = rShoulder * Math.sin(theta);

    // frustrum
    // const frustumEdgeX2 = noseH + (rBody - noseL - shoulderH) / Math.tan(theta); // X-coordinate where frustum edge meets shoulder [m] //UNUSED
    const frustumEdgeR2 = rBody - shoulderH; // R-coordinate where frustum edge meets shoulder [m]
    const frustumSide = Math.sqrt((rBody - noseL) * (rBody - noseL) + (frustumEdgeR2 - noseL) * (frustumEdgeR2 - noseL)); // Q: Is this correct?
    // const frustumArea = Math.PI * (rBody + noseL) * frustumSide; // UNUSED
    const frustumL = rBody - noseL - shoulderH; // length
    const frustumH = frustumL / Math.tan(theta); // height
    // const gradient = frustumH / frustumL; //UNUSED
    // shoulderStartR and shoulderStartX same as frustumEdgeX2 and frustumEdgeR2
    // cone (if the geometry is extended)
    const coneH = (rBody - shoulderH) / Math.tan(theta); // cone height
    const coneVertex = -coneH + noseH + frustumH;

    // define total number of mesh points for nose, frustum and shoulder
    const { nop } = this.mesh;
    // define number of mesh divisions depending on the lengths of the three parts
    const curveLength = rNose * (Math.PI / 2 - theta) + frustumSide + rShoulder * theta;
    const nopNose = Math.round(nop * rNose * (Math.PI / 2 - theta) / curveLength);
    const nopFrustum = Math.round(nop * frustumSide / curveLength);
    const nopShoulder = nop - nopNose - nopFrustum + 2; // two mesh points are doubled later, add to shoulder
    // define radial vectors for nose (r1), frustum (r2) and shoulder (r3)
    const r1 = numeric.linspace(0, noseL, nopNose);
    const r2 = numeric.linspace(noseL, frustumL + noseL, nopFrustum).slice(1);
    const r3 = numeric.linspace(frustumL + noseL, rBody, nopShoulder).slice(1);

    // define rhoC and thetaC coordinates (cylindrical coordinates)
    const rhoC = r1.concat(r2, r3); // combine to create single vector for full heatshield radius - [r1,r2(2:end),r3(2:end)]
    const thetaC = numeric.linspace(0, 2 * Math.PI, rhoC.length); // theta angle vector for full heatshield

    // define rLath, xLath coordinates for LatheGeometry in THREE.js
    const rLath = rhoC; // r
    const xLath = r1.map(i => -Math.sqrt(Math.abs(rNose * rNose - i * i)) + rNose).concat(
      r2.map(i => i / Math.tan(theta) + coneVertex),
      r3.map(i => -Math.sqrt(Math.abs(rShoulder * rShoulder - (i - (rBody - rShoulder)) * (i - (rBody - rShoulder)))) + (shoulderL + frustumH + noseH)),
    ); // x

    // call meshgrid function for theta and rho
    const [thetaCM, rhoCM] = this.meshgrid(thetaC, rhoC);

    // define x matrix
    const x = [];
    for (let i = 0; i < rhoCM.length; i++) {
      const temp = [];
      for (let j = 0; j < rhoCM.length; j++) {
        if (rhoCM[i][j] <= noseL) {
          temp.push(-Math.sqrt(Math.abs(rNose * rNose - rhoCM[i][j] * rhoCM[i][j])) + rNose);
        } else if (rhoCM[i][j] > noseL && rhoCM[i][j] <= (frustumL + noseL)) {
          temp.push(rhoCM[i][j] / Math.tan(theta) + coneVertex);
        } else if (rhoCM[i][j] <= rBody) {
          temp.push(-Math.sqrt(Math.abs(rShoulder * rShoulder - (rhoCM[i][j] - (rBody - rShoulder)) * (rhoCM[i][j] - (rBody - rShoulder)))) + (shoulderL + frustumH + noseH));
        }
      }
      // temp = rhoCM[i].map(function(item) {return -Math.sqrt(Math.abs(rNose*rNose-item*itemp))+rNose;});
      x.push(temp);
    }

    const xN = x.map(i => i.map(j => -j));

    // y and z matrices with cartesian coordinates of the meshgrid
    const y = numeric.mul(rhoCM, numeric.cos(thetaCM)); // y = rhoC.*cos(thetaC)
    const z = numeric.mul(rhoCM, numeric.sin(thetaCM)); // z = rhoC.*sin(thetaC)

    // set the results into lander.mesh
    this.mesh.x = xN;
    this.mesh.y = y;
    this.mesh.z = z;
    this.mesh.rLath = rLath;
    this.mesh.xLath = xLath;
  },
};
lander.defineMesh();
// /////////////////
// TRAJ (MODULE 2)//
// /////////////////

// define Mars constants (global object)
const Mars = {
  r: 3396190, // Mars equatorial radius [m]
  omega: 7.08824 * (10 ** -5), // Mars rotation rate [rad/s]
  mu: 4.282837 * (10 ** 13), // Mars standard gravitational parameter - GM - [m^3/s^2]
  J2: 1.960454 * (10 ** -3), // Mars second dynamic factor
  // coefficients to calculate convective and radiative heat rate
  k: 1.898 * (10 ** -4), // Mars Sutton&Graves coeff.
  C: 2.35 * (10 ** 4), // Mars Tauber&Sutton coeff. C
  a: 0.526, // Mars Tauber&Sutton coeff. a
  b: 1.19, // Mars Tauber&Sutton coeff. b
};
Object.freeze(Mars); // shallow freeze the object

// define trajectory options module that can be changed by the gui
const trajectoryOptions = {
  // define initial conditions
  cond: {
    // initial angles
    FPADeg: 17, // * Flight Path Angle [deg down from horizontal]
    aziDeg: 77, // * azimuth [deg from N]
    latDeg: 69.36, // * latitude [deg N]
    lonDeg: 197.69, // * longitude [deg E]
    alphaDeg: -11.1, // * angle of attack [deg]
    betaDeg: 0, // * angle of sideslip [deg]
    bankDeg: 0, // bank angle [deg]

    // convert all angles to radians
    get FPA() { return this.FPADeg * Math.PI / 180; }, // * Flight Path Angle [rad]
    get azi() { return this.aziDeg * Math.PI / 180; }, // * azimuth [rad]
    get lat() { return this.latDeg * Math.PI / 180; }, // * latitude [rad]
    get lon() { return this.lonDeg * Math.PI / 180; }, // * longitude [rad]
    get alpha() { return this.alphaDeg * Math.PI / 180; }, // * angle of attack [rad]
    get beta() { return this.betaDeg * Math.PI / 180; }, // * angle of sideslip [rad]
    get bank() { return this.bankDeg * Math.PI / 180; }, // bank angle [rad]

    // define other initial angles
    get alphaTotal() { return Math.acos(Math.cos(this.alpha) * Math.cos(this.beta)); }, // total angle of attack [rad]
    get phi_a() { return Math.atan2(Math.tan(this.beta), Math.sin(this.alpha)); }, // aerodynamic roll angle [rad]

    // initial roll, pitch, yaw
    get roll() { return Math.atan2(Math.sin(this.bank), (Math.cos(this.alpha) * (Math.cos(this.bank) - Math.tan(-this.FPA) * Math.tan(this.alpha)))); }, // always 0 for bank = 0, roll angle [rad] (Wagner, 1965)
    get pitch() { return this.alpha - this.FPA; }, // * pitch angle [rad] (FPA down is positive)
    get yaw() { return this.azi - this.beta; }, // * yaw angle [rad]

    // initial rotational speed
    omegaXDeg: 0.05, // roll rate, p [deg/s]
    omegaYDeg: 0.2, // pitch rate, q [deg/s]
    omegaZDeg: 0.05, // yaw rate, r [deg/s]
    // convert to rad/s
    get omegaX() { return this.omegaXDeg * Math.PI / 180; }, // roll rate, p [rad/s]
    get omegaY() { return this.omegaYDeg * Math.PI / 180; }, // pitch rate, q [rad/s]
    get omegaZ() { return this.omegaZDeg * Math.PI / 180; }, // yaw rate, r [rad/s]
    get omega() { return [this.omegaX, this.omegaY, this.omegaZ]; }, // initial angular velocity vector [rad/s]

    // define initial position
    h: 125000, // altitude [m]
    get r() { return Mars.r + this.h; }, // distance between the lander and Mars centre [m]

    // define initial velocity
    vInit: 7400, // * relative entry velocity [m/s]
    vInf: 7400, // * relative wind velocity [m/s]

    get u() { return this.vInf * Math.cos(this.FPA) * Math.cos(this.azi); }, // * component of velocity in N direction [m/s]
    get v() { return this.vInf * Math.cos(this.FPA) * Math.sin(this.azi); }, // * component of velocity in E direction [m/s]
    get w() { return this.vInf * Math.sin(this.FPA); }, // * component of velocity in D direction [m/s]

    // define initial attitude (quaternions)
    get e0() { return Math.cos(this.roll / 2) * Math.cos(this.pitch / 2) * Math.cos(this.yaw / 2) + Math.sin(this.roll / 2) * Math.sin(this.pitch / 2) * Math.sin(this.yaw / 2); },
    get e1() { return Math.sin(this.roll / 2) * Math.cos(this.pitch / 2) * Math.cos(this.yaw / 2) - Math.cos(this.roll / 2) * Math.sin(this.pitch / 2) * Math.sin(this.yaw / 2); },
    get e2() { return Math.cos(this.roll / 2) * Math.sin(this.pitch / 2) * Math.cos(this.yaw / 2) + Math.sin(this.roll / 2) * Math.cos(this.pitch / 2) * Math.sin(this.yaw / 2); },
    get e3() { return Math.cos(this.roll / 2) * Math.cos(this.pitch / 2) * Math.sin(this.yaw / 2) - Math.sin(this.roll / 2) * Math.sin(this.pitch / 2) * Math.cos(this.yaw / 2); },

    // other unused variables
    get vOrb() { return Math.sqrt(Mars.mu / (Mars.r + this.h)); }, // orbital velocity for initial altitude
    get tOrb() { return 2 * Math.PI * (Mars.r + this.h) / this.vOrb; }, // orbital period for circular orbit
    Qload: 0, // heat load at entry interface point (EIP)

  },
  // define mission object containing objects with atmospheric properties (object)
  // store current
  current: 'pathfinder',
  // define TauberSuttonVelocityFunction (object)
  TauberSuttonVelocityFunction: {
    velocity: [6000, 6150, 6300, 6500, 6700, 6900, 7000, 7200, 7400, 7600, 7800, 8000, 8200, 8400, 8600, 8800, 9000],
    fv: [0.20, 1, 1.950, 3.420, 5.100, 7.100, 8.100, 10.20, 12.50, 14.80, 17.10, 19.20, 21.40, 24.10, 26, 28.90, 32.80],
  },
  // define mission object containing objects with atmospheric properties (object)
  mission: {
    // Pathfinder atmosphere
    pathfinder: {
      altitude: [0, 7352.94, 14705.9, 22058.8, 29411.8, 36764.7, 44117.6, 51470.6, 58823.5, 66176.5, 73529.4, 80882.4, 88235.3, 95588.2, 102941, 110294, 117647, 125000, 132353, 139706, 147059, 154412, 161765, 169118, 176471, 183824, 191176, 198529, 205882, 213235, 220588, 227941, 235294, 242647, 250000],
      temperature: [199.85, 213.24, 198.62, 186.28, 178.80, 172.26, 162.96, 154.80, 148.80, 139.93, 135.03, 136.55, 127.42, 122.04, 119.14, 115.16, 119.08, 129.14, 138.40, 144.20, 147.42, 149.12, 150.02, 150.47, 150.70, 150.85, 150.95, 150.99, 151.03, 151.04, 151.06, 151.07, 151.08, 151.08, 151.08],
      pressure: [687.85, 357.27, 177.99, 84.724, 38.856, 17.371, 7.4808, 3.0919, 1.2279, 0.46788, 0.16976, 0.061781, 0.021902, 0.0073144, 0.0024159, 0.00079147, 0.00026671, 9.8606e-05, 4.0940e-05, 1.8607e-05, 9.0356e-06, 4.6453e-06, 2.5136e-06, 1.4316e-06, 8.5982e-07, 5.3453e-07, 3.4166e-07, 2.2868e-07, 1.5468e-07, 1.0822e-07, 7.6339e-08, 5.5418e-08, 4.0282e-08, 3.0242e-08, 2.2704e-08],
      density: [0.017639, 0.0087843, 0.0046950, 0.0023827, 0.0011385, 0.00052866, 0.00024033, 0.00010494, 4.3278e-05, 1.7540e-05, 6.6336e-06, 2.3732e-06, 8.9967e-07, 3.1204e-07, 1.0412e-07, 3.4629e-08, 1.1087e-08, 3.6569e-09, 1.3616e-09, 5.6871e-10, 2.5638e-10, 1.2264e-10, 6.1546e-11, 3.2504e-11, 1.8220e-11, 1.0515e-11, 6.2124e-12, 3.9405e-12, 2.4894e-12, 1.6646e-12, 1.1086e-12, 7.7363e-13, 5.3193e-13, 3.8653e-13, 2.7736e-13],
      cp: [746.22, 762.84, 744.66, 728.61, 718.55, 709.51, 696.28, 684.26, 675.17, 661.29, 653.46, 656.31, 642.21, 635.86, 636.59, 639.27, 659.90, 693.32, 727.94, 761.28, 795.93, 833.99, 876.34, 922.45, 969.18, 1016.23, 1063.78, 1102.08, 1141.22, 1170.08, 1198.81, 1218.12, 1238.98, 1250.23, 1262.54],
      viscosity: [1.0237e-05, 1.0945e-05, 1.0171e-05, 9.5105e-06, 9.1057e-06, 8.7498e-06, 8.2413e-06, 7.7915e-06, 7.4595e-06, 6.9660e-06, 6.6931e-06, 6.7814e-06, 6.2846e-06, 6.0284e-06, 5.9642e-06, 5.8825e-06, 6.2676e-06, 7.0180e-06, 7.7368e-06, 8.2766e-06, 8.6926e-06, 9.0454e-06, 9.3724e-06, 9.6859e-06, 9.9781e-06, 1.0254e-05, 1.0516e-05, 1.0721e-05, 1.0921e-05, 1.1068e-05, 1.1210e-05, 1.1307e-05, 1.1408e-05, 1.1465e-05, 1.1525e-05],
    },
    // Curiosity atmosphere
    curiosity: {
      altitude: [0, 7352.94, 14705.9, 22058.8, 29411.8, 36764.7, 44117.6, 51470.6, 58823.5, 66176.5, 73529.4, 80882.4, 88235.3, 95588.2, 102941, 110294, 117647, 125000, 132353, 139706, 147059, 154412, 161765, 169118, 176471, 183824, 191176, 198529, 205882, 213235, 220588, 227941, 235294, 242647, 250000],
      temperature: [194.541, 218.893, 200.782, 180.643, 173.902, 175.245, 174.950, 169.970, 162.451, 150.246, 140.370, 141.065, 145.247, 148.469, 142.646, 131.164, 125.699, 126.025, 131.041, 137.403, 142.069, 144.821, 146.278, 147.021, 147.406, 147.616, 147.725, 147.785, 147.828, 147.850, 147.867, 147.876, 147.884, 147.890, 147.895],
      pressure: [751.267, 392.079, 198.590, 93.8924, 41.9237, 18.6125, 8.31677, 3.68289, 1.59161, 0.650339, 0.249292, 0.0924766, 0.0351107, 0.0138635, 0.00549805, 0.00206911, 0.000745505, 0.000267921, 0.000100236, 4.02404e-05, 1.73348e-05, 7.91875e-06, 3.82610e-06, 1.95609e-06, 1.05771e-06, 6.02104e-07, 3.61405e-07, 2.28891e-07, 1.48725e-07, 1.00708e-07, 6.98071e-08, 4.96649e-08, 3.59762e-08, 2.67246e-08, 1.98558e-08],
      density: [0.0195581, 0.00938476, 0.00518105, 0.00272454, 0.00126500, 0.000558136, 0.000249716, 0.000113644, 5.14132e-05, 2.26365e-05, 9.31522e-06, 3.47447e-06, 1.28318e-06, 4.88830e-07, 1.99619e-07, 8.07924e-08, 3.01452e-08, 1.06699e-08, 3.77249e-09, 1.41003e-09, 5.69444e-10, 2.44570e-10, 1.10869e-10, 5.28509e-11, 2.64369e-11, 1.37981e-11, 7.60120e-12, 4.46644e-12, 2.66222e-12, 1.68757e-12, 1.09031e-12, 7.30914e-13, 4.96918e-13, 3.51821e-13, 2.45314e-13],
      cp: [739.399, 769.579, 747.362, 721.057, 711.818, 713.673, 713.265, 706.319, 695.579, 677.442, 662.113, 663.669, 671.423, 678.380, 672.411, 657.959, 655.080, 664.148, 682.180, 705.200, 730.312, 758.366, 791.427, 830.816, 876.441, 927.381, 981.024, 1031, 1081.27, 1122.73, 1159.83, 1189.21, 1214.72, 1231.61, 1249.41],
      viscosity: [9.96245e-06, 1.12481e-05, 1.02939e-05, 9.21212e-06, 8.84420e-06, 8.91658e-06, 8.90014e-06, 8.62836e-06, 8.21605e-06, 7.54237e-06, 6.99393e-06, 7.03817e-06, 7.29187e-06, 7.51471e-06, 7.25498e-06, 6.67605e-06, 6.45446e-06, 6.56705e-06, 6.95356e-06, 7.44290e-06, 7.87025e-06, 8.22612e-06, 8.54938e-06, 8.87181e-06, 9.20401e-06, 9.54308e-06, 9.87378e-06, 1.01642e-05, 1.04386e-05, 1.06569e-05, 1.08445e-05, 1.09911e-05, 1.11152e-05, 1.11987e-05, 1.12842e-05],
    },
    // Opportunity atmosphere
    opportunity: {
      altitude: [0, 7352.94, 14705.9, 22058.8, 29411.8, 36764.7, 44117.6, 51470.6, 58823.5, 66176.5, 73529.4, 80882.4, 88235.3, 95588.2, 102941, 110294, 117647, 125000, 132353, 139706, 147059, 154412, 161765, 169118, 176471, 183824, 191176, 198529, 205882, 213235, 220588, 227941, 235294, 242647, 250000],
      temperature: [199.76, 213.32, 201.18, 188.59, 180.64, 179.74, 179.18, 174.72, 165.02, 152.48, 139.03, 133.70, 133.44, 136.99, 130.46, 115.72, 115.36, 131.56, 145.96, 154.77, 159.78, 162.50, 163.97, 164.78, 165.24, 165.51, 165.66, 165.75, 165.82, 165.84, 165.87, 165.89, 165.90, 165.91, 165.91],
      pressure: [647.761, 338.553, 169.540, 81.5839, 37.7004, 17.1878, 7.8465, 3.5570, 1.5663, 0.64930, 0.25020, 0.090600, 0.032100, 0.011600, 0.0042, 0.0014, 0.00044401, 0.00015369, 6.1701e-05, 2.7138e-05, 1.2710e-05, 6.2590e-06, 3.2298e-06, 1.7482e-06, 9.9323e-07, 5.9002e-07, 3.6837e-07, 2.4330e-07, 1.6356e-07, 1.1791e-07, 8.4995e-08, 6.5429e-08, 5.0775e-08, 3.9781e-08, 3.3218e-08],
      density: [0.0165, 0.00833, 0.00442, 0.00227, 0.00110, 0.000502, 0.000230, 0.000107, 4.99e-05, 2.23e-05, 9.43e-06, 3.57e-06, 1.27e-06, 4.46e-07, 1.68e-07, 6.39e-08, 1.99e-08, 5.98e-09, 2.12e-09, 8.55e-10, 3.76e-10, 1.75e-10, 8.54e-11, 4.35e-11, 2.30e-11, 1.26e-11, 7.21e-12, 4.39e-12, 2.65e-12, 1.78e-12, 1.16e-12, 8.30e-13, 5.91e-13, 4.15e-13, 3.27e-13],
      cp: [746.15, 763.04, 747.96, 731.69, 721.04, 719.81, 719.05, 712.91, 699.19, 680.67, 659.69, 651.08, 650.93, 657.50, 648.27, 625.02, 629.83, 664.95, 696.35, 720.37, 741.65, 763.80, 789.38, 819.86, 855.51, 896.19, 941.07, 984.94, 1034.05, 1070.69, 1114.73, 1140.88, 1168.59, 1198.88, 1209.71],
      viscosity: [1.0221e-05, 1.0937e-05, 1.0296e-05, 9.6228e-06, 9.1935e-06, 9.1437e-06, 9.1131e-06, 8.8715e-06, 8.3422e-06, 7.6524e-06, 6.9060e-06, 6.6093e-06, 6.5988e-06, 6.8094e-06, 6.4753e-06, 5.6857e-06, 5.7560e-06, 6.7782e-06, 7.6956e-06, 8.3099e-06, 8.7323e-06, 9.0569e-06, 9.3472e-06, 9.6377e-06, 9.9399e-06, 1.0254e-05, 1.0573e-05, 1.0863e-05, 1.1164e-05, 1.1377e-05, 1.1618e-05, 1.1756e-05, 1.1898e-05, 1.2046e-05, 1.2099e-05],
    },
    // Phoenix atmosphere
    phoenix: {
      altitude: [0, 7352.94, 14705.9, 22058.8, 29411.8, 36764.7, 44117.6, 51470.6, 58823.5, 66176.5, 73529.4, 80882.4, 88235.3, 95588.2, 102941, 110294, 117647, 125000, 132353, 139706, 147059, 154412, 161765, 169118, 176471, 183824, 191176, 198529, 205882, 213235, 220588, 227941, 235294, 242647, 250000],
      temperature: [191.864, 209.266, 192.737, 175.560, 161.970, 155.099, 153.433, 155.459, 159.757, 161.637, 158.347, 151.329, 141.997, 136.794, 135.898, 128.260, 120.140, 128.648, 152.137, 174.895, 192.620, 204.952, 212.408, 217.354, 220.621, 222.317, 223.415, 224.157, 224.519, 224.835, 224.980, 225.114, 225.192, 225.248, 225.298],
      pressure: [845.146, 433.364, 212.077, 97.5030, 41.8197, 17.0780, 6.83704, 2.75275, 1.13534, 0.478, 0.201470, 0.0826777, 0.0325379, 0.0123104, 0.00462, 0.0017, 0.00059, 0.000208, 8.42e-05, 3.98e-05, 2.06e-05, 1.14e-05, 6.67e-06, 4.02e-06, 2.49e-06, 1.59e-06, 1.05e-06, 6.99e-07, 4.86e-07, 3.39e-07, 2.47e-07, 1.81e-07, 1.36e-07, 1.03e-07, 7.87e-08],
      density: [0.0223, 0.0109, 0.00578, 0.00291, 0.00135, 0.000578, 0.000234, 9.30e-05, 3.74e-05, 1.55e-05, 6.67e-06, 2.85e-06, 1.19e-06, 4.67e-07, 1.75e-07, 6.76e-08, 2.49e-08, 8.11e-09, 2.74e-09, 1.11e-09, 5.10e-10, 2.60e-10, 1.43e-10, 8.16e-11, 4.78e-11, 2.92e-11, 1.82e-11, 1.15e-11, 7.62e-12, 4.99e-12, 3.47e-12, 2.39e-12, 1.71e-12, 1.24e-12, 8.94e-13],
      cp: [736.091, 758.315, 737.191, 714.048, 694.720, 684.582, 682.121, 685.238, 691.709, 694.597, 690.009, 680.077, 667.130, 661.334, 662.217, 652.718, 644.558, 665.395, 706.485, 743.133, 772.102, 795.126, 814.570, 833.321, 852.489, 873.287, 896.111, 921.389, 948.247, 978.036, 1006.57, 1038.09, 1066.46, 1094.90, 1124.25],
      viscosity: [9.76e-06, 1.07e-05, 9.82e-06, 8.91e-06, 8.17e-06, 7.79e-06, 7.70e-06, 7.82e-06, 8.05e-06, 8.16e-06, 7.98e-06, 7.60e-06, 7.12e-06, 6.87e-06, 6.86e-06, 6.48e-06, 6.12e-06, 6.71e-06, 8.13e-06, 9.47e-06, 1.05e-05, 1.13e-05, 1.18e-05, 1.22e-05, 1.26e-05, 1.28e-05, 1.31e-05, 1.33e-05, 1.36e-05, 1.38e-05, 1.40e-05, 1.42e-05, 1.44e-05, 1.46e-05, 1.48e-05],
    },
    // Schiaparelli atmosphere
    schiaparelli: {
      altitude: [0, 7352.94, 14705.9, 22058.8, 29411.8, 36764.7, 44117.6, 51470.6, 58823.5, 66176.5, 73529.4, 80882.4, 88235.3, 95588.2, 102941, 110294, 117647, 125000, 132353, 139706, 147059, 154412, 161765, 169118, 176471, 183824, 191176, 198529, 205882, 213235, 220588, 227941, 235294, 242647, 250000],
      temperature: [284.317, 226.882, 218.725, 209.186, 194.974, 179.969, 173.044, 165.429, 157.696, 160.376, 161.386, 151.818, 141.966, 136.419, 133.411, 124.331, 105.328, 107.005, 157.355, 204.013, 236.098, 254.979, 268.309, 277.695, 282.318, 286.052, 288.199, 289.647, 290.619, 291.216, 291.654, 291.938, 292.136, 292.294, 292.373],
      pressure: [725.89, 392.29, 206.19, 106, 52.533, 24.586, 11.058, 4.8234, 2.0243, 0.83980, 0.35521, 0.14693, 0.057706, 0.021596, 0.0078835, 0.0027719, 0.00085874, 0.00023410, 8.4119e-05, 4.1061e-05, 2.2846e-05, 1.3766e-05, 8.5551e-06, 5.4340e-06, 3.5607e-06, 2.3551e-06, 1.5839e-06, 1.0822e-06, 7.4745e-07, 5.2743e-07, 3.7557e-07, 2.7365e-07, 2.0197e-07, 1.5109e-07, 1.1608e-07],
      density: [0.015268, 0.0090745, 0.0049472, 0.0026574, 0.0014136, 0.00071569, 0.00033534, 0.00015276, 6.7365e-05, 2.7570e-05, 1.1546e-05, 5.0649e-06, 2.1284e-06, 8.3041e-07, 3.0946e-07, 1.1662e-07, 4.2449e-08, 1.1568e-08, 2.7763e-09, 1.0333e-09, 4.9088e-10, 2.7348e-10, 1.5907e-10, 9.5925e-11, 6.1176e-11, 3.9171e-11, 2.5631e-11, 1.7014e-11, 1.1377e-11, 7.7734e-12, 5.3239e-12, 3.7314e-12, 2.6375e-12, 1.8775e-12, 1.3838e-12],
      cp: [841.70, 779.28, 769.63, 758.02, 740.06, 720.12, 710.55, 699.74, 688.42, 692.41, 693.97, 679.80, 664.62, 655.88, 651.20, 636.22, 602.70, 610.75, 696.91, 760.83, 800.15, 823.03, 839.73, 852.41, 861.54, 870.39, 878.93, 888.21, 898.28, 910, 923.01, 937.81, 954.39, 972.32, 991.59],
      viscosity: [1.4521e-05, 1.1639e-05, 1.1214e-05, 1.0713e-05, 9.9581e-06, 9.1500e-06, 8.7741e-06, 8.3585e-06, 7.9341e-06, 8.0828e-06, 8.1404e-06, 7.6148e-06, 7.0701e-06, 6.7636e-06, 6.5983e-06, 6.0926e-06, 5.0387e-06, 5.2168e-06, 8.1228e-06, 1.0687e-05, 1.2390e-05, 1.3387e-05, 1.4098e-05, 1.4614e-05, 1.4910e-05, 1.5171e-05, 1.5369e-05, 1.5548e-05, 1.5717e-05, 1.5886e-05, 1.6060e-05, 1.6239e-05, 1.6426e-05, 1.6616e-05, 1.6804e-05],
    },
    // Spirit atmosphere
    spirit: {
      altitude: [0, 7352.94, 14705.9, 22058.8, 29411.8, 36764.7, 44117.6, 51470.6, 58823.5, 66176.5, 73529.4, 80882.4, 88235.3, 95588.2, 102941, 110294, 117647, 125000, 132353, 139706, 147059, 154412, 161765, 169118, 176471, 183824, 191176, 198529, 205882, 213235, 220588, 227941, 235294, 242647, 250000],
      temperature: [189.37, 219.83, 205.63, 193.52, 182.42, 175.80, 170.71, 165.91, 161.49, 155.52, 141.59, 134.58, 131.23, 127.19, 120.91, 110.72, 109.90, 126.91, 142.74, 151.55, 156.27, 158.72, 159.97, 160.63, 160.99, 161.17, 161.27, 161.35, 161.38, 161.41, 161.42, 161.44, 161.44, 161.45, 161.45],
      pressure: [685.01, 358.53, 182.8, 89.39, 41.96, 19.03, 8.43, 3.66, 1.56, 0.650, 0.250, 0.0933, 0.0329, 0.0114, 0.00383, 0.00119, 0.000351, 0.000116, 4.43e-05, 1.90e-05, 8.71e-06, 4.19e-06, 2.11e-06, 1.12e-06, 6.23e-07, 3.68e-07, 2.30e-07, 1.49e-07, 1.02e-07, 7.20e-08, 5.31e-08, 3.95e-08, 3.10e-08, 2.43e-08, 1.90e-08],
      density: [0.0181, 0.00856, 0.00466, 0.00242, 0.00121, 0.000568, 0.000259, 0.000116, 5.08e-05, 2.18e-05, 9.42e-06, 3.66e-06, 1.32e-06, 4.72e-07, 1.65e-07, 5.62e-08, 1.66e-08, 4.73e-09, 1.56e-09, 6.18e-10, 2.66e-10, 1.21e-10, 5.78e-11, 2.87e-11, 1.47e-11, 8.00e-12, 4.61e-12, 2.70e-12, 1.69e-12, 1.09e-12, 7.38e-13, 4.91e-13, 3.64e-13, 2.64e-13, 1.92e-13],
      cp: [732.72, 770.92, 753.60, 738.16, 723.45, 714.40, 707.29, 700.45, 694.04, 685.22, 663.74, 652.47, 647.19, 640.96, 631.22, 615.03, 618.81, 655.66, 689.12, 713.36, 734.75, 757.41, 784.10, 816.36, 854.88, 899.20, 945.67, 995.34, 1041.33, 1085.42, 1121.07, 1160.91, 1181.37, 1205.31, 1226.07],
      viscosity: [9.66e-06, 1.13e-05, 1.05e-05, 9.88e-06, 9.29e-06, 8.93e-06, 8.65e-06, 8.39e-06, 8.15e-06, 7.82e-06, 7.05e-06, 6.66e-06, 6.47e-06, 6.25e-06, 5.92e-06, 5.38e-06, 5.43e-06, 6.50e-06, 7.50e-06, 8.11e-06, 8.52e-06, 8.83e-06, 9.12e-06, 9.41e-06, 9.73e-06, 1.01e-05, 1.04e-05, 1.07e-05, 1.10e-05, 1.12e-05, 1.14e-05, 1.16e-05, 1.17e-05, 1.18e-05, 1.19e-05],
    },
  },
  /**
   * Calculates atmospheric properties from known atmospheric data.
   *
   * @param {number} h - current altitude of the lander
   * @param {number} vInf - current lander velocity
   * @param {object} lander - contains dimensions of the lander
   * @param {object} mission - contains atmospheric data
   *
   * @return {array} - array with atmospheric properties
   */
  calcAtmoProperties(h, vInf, lander) {
    // variables from lander object
    const { diameter } = lander;

    // atmosphere conditions
    const atm = this.mission[trajectoryOptions.current];

    let tInf; let pInf; let densInf; let cpInf; let viscInf;

    if (h > 250000) {
      // assuming the length of the data is 35
      tInf = atm.temperature[34]; // freestream temperature [K]
      pInf = atm.pressure[34]; // freestream pressure [Pa]
      densInf = atm.density[34]; // freestream density [kg/m^3]
      cpInf = atm.cp[34]; // specific heat at constant pressure [J/kgK]
      viscInf = atm.viscosity[34]; // dynamic viscosity [Ns/m2]
    } else {
      tInf = numeric.spline(atm.altitude, atm.temperature).at(h);
      pInf = numeric.spline(atm.altitude, atm.pressure).at(h);
      densInf = numeric.spline(atm.altitude, atm.density).at(h);
      cpInf = numeric.spline(atm.altitude, atm.cp).at(h);
      viscInf = numeric.spline(atm.altitude, atm.viscosity).at(h);
    }

    // calculate freestream properties
    const R = pInf / (densInf * tInf); // gas constant R
    const cvInf = cpInf - R; // specific heat at constant volume Cv
    const gamma = cpInf / cvInf; // ratio of specific heats
    const aSound = Math.sqrt(gamma * R * tInf); // speed of sound
    const qInf = 0.5 * densInf * vInf * vInf; // dynamic pressure
    const mach = vInf / aSound; // Mach number
    const Re = densInf * vInf * diameter / viscInf; // Reynolds Number
    const Kn = mach / Re * Math.sqrt(Math.PI * gamma / 2); // Knudsen Number
    const cpMax = (2 / (mach * mach * gamma)) * (((((gamma + 1) * (gamma + 1) * mach * mach) / (4 * gamma * mach * mach - 2 * (gamma - 1))) ** (gamma / (gamma - 1))) * ((1 - gamma + 2 * gamma * mach * mach) / (gamma + 1)) - 1);
    // return atmospheric properties
    return [tInf, pInf, densInf, cpInf, viscInf, R, cvInf, gamma, aSound, qInf, mach, Re, Kn, cpMax];
  },
};

// ////////////////////
// RESULTS (MODULE 3)//
// ////////////////////
// define a worker that will calculate new trajectory in a new thread
const worker = new Worker('scripts/trajectory-worker.js');
// add event listener to the worker on the main thread
worker.addEventListener('message', (event) => {
  // register results sent by the worker
  results.sol = event.data.data1;

  // transform the raw results into arrays used by the animation system (THREE.js)
  const { scale } = lander.mesh;
  results.createArrays(results.sol, results.arr, scale);
  // create spline objects for interpolation
  results.createSplines(results.arr, results.spline);

  // initialize live outputs and animation for three.js
  liveOutputs.init();
  G.defineAnim(results);

  // hide please wait notification
  const notif = document.querySelector('#pleasewait');
  notif.style.visibility = 'hidden';

  // display notification to notify the user that the trajectory has been calculated
  const trajCalculated = document.querySelector('#trajCalculated');
  trajCalculated.style.visibility = 'visible';
});

// This module transforms results obtained from traj module into
// a format that THREE.js can use. Create splines that are used for interpolation.
const results = {
  sol: null, // raw solution from trajectory module
  arr: {}, // arrays used by the animation system (THREE.js)
  spline: {}, // spline objects for interpolation
  /**
  * Container function for trajectorySolve.solve(), createArrays and createSplines functions.
  * This runs when the button "Calculate" is pressed
  */
  run() {
    // show "please wait" notification
    const notif = document.querySelector('#pleasewait');
    notif.style.visibility = 'visible';

    // define data that will be sent to a worker
    // trajectory data
    const data1 = {
      cond: trajectoryOptions.cond,
      current: trajectoryOptions.current,
    };
    // lander geometry data
    const data2 = {
      mass: lander.mass,
      diameter: lander.diameter,
      rBody: lander.rBody,
      rNose: lander.rNose,
      rShoulder: lander.rShoulder,
      sphereConeAngleDeg: lander.sphereConeAngleDeg,
      theta: lander.theta,
      refArea: lander.refArea,
      COMx: lander.COMx,
      COMpayload: lander.COMpayload,
      rPayload: lander.rPayload,
      heightPayload: lander.heightPayload,
      hsThickness: lander.hsThickness,
      hsDensity: lander.hsDensity,
      mesh: lander.mesh,
    };
    // send a message with all necessary data to a worker to calculate trajectory
    worker.postMessage({
      data1,
      data2,
    });
  },
  /**
  * Transforms initial sol object from ODE solver into separate arrays.
  *
  * @param {object} sol - initial results from ODE solver
  * @param {object} arr - empty object for the results
  * @param {number} scale - scales all results for graphics
  *
  */
  createArrays(sol, arr, scale) {
    // initialize
    const quaternionNED = new THREE.Quaternion();
    const rotationNED = new THREE.Matrix4();
    const solLength = sol.x.length;

    // preallocate arrays
    const time = sol.x;
    const alt = new Array(solLength);
    const pos = new Array(3 * solLength);
    const rotN = new Array(4 * solLength);
    const rotB = new Array(4 * solLength);
    const xVel = new Array(solLength);
    const yVel = new Array(solLength);
    const zVel = new Array(solLength);
    const vel = new Array(solLength);
    const qInf = new Array(solLength);
    const tempNose = new Array(solLength);
    const AoA = new Array(solLength);

    // create arrays
    for (let i = 0; i < solLength; i++) {
      const [r, lat, lon, u, v, w, e0, e1, e2, e3] = sol.y[i];


      // altitude (live outputs) not scaled
      alt[i] = r - Mars.r;

      // length of velocity vector (live outputs, graphics) not scaled
      vel[i] = Math.sqrt(u * u + v * v + w * w);

      // dynamic pressure (live outputs)
      qInf[i] = trajectoryOptions.calcAtmoProperties(alt[i], vel[i], lander)[9];

      // angle of attack
      const vInfGc = [u, v, w];
      const Ggd2b = [[e0 * e0 + e1 * e1 - e2 * e2 - e3 * e3, 2 * (e1 * e2 + e0 * e3), 2 * (e1 * e3 - e0 * e2)], [2 * (e1 * e2 - e0 * e3), e0 * e0 - e1 * e1 + e2 * e2 - e3 * e3, 2 * (e0 * e1 + e2 * e3)], [2 * (e1 * e3 + e0 * e2), 2 * (e2 * e3 - e0 * e1), e0 * e0 - e1 * e1 - e2 * e2 + e3 * e3]];
      const Ggc2gd = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
      const Ggc2b = numeric.dot(Ggc2gd, Ggd2b);
      const vInfB = numeric.dot(Ggc2b, vInfGc); // velocity vector in body frame MxV
      const [vX, , vZ] = vInfB; // x, y, z body frame velocity coordinates
      AoA[i] = Math.atan2(vZ, vX); // angle of attack

      // calculate temperature at nose (live outputs)
      const densInf = trajectoryOptions.calcAtmoProperties(alt[i], vel[i], lander)[2];
      const qConv = Mars.k * Math.sqrt(densInf / lander.rNose) * vel[i] * vel[i] * vel[i] / 10000; // Sutton-Graves relation for convective heat flux at stagnation point [W/cm2]
      let qRad;
      if (vel[i] > 6500) {
        const TSv = numeric.spline(trajectoryOptions.TauberSuttonVelocityFunction.velocity, trajectoryOptions.TauberSuttonVelocityFunction.fv).at(vel[i]);
        qRad = Mars.C * (lander.rNose ** Mars.a) * (densInf ** Mars.b) * TSv; // Tauber-Sutton relation for radiative heat flux at stagnation point [W/cm2]
      } else {
        qRad = 0;
      }
      const qTotal = qConv + qRad;
      tempNose[i] = (qTotal / (5.670373e-12)) ** 0.25;

      // direction of velocity vector (graphics) not scaled
      xVel[i] = u / vel[i]; // x
      yVel[i] = v / vel[i]; // y
      zVel[i] = w / vel[i]; // z

      // position NED, (graphics) scaled
      pos[3 * i] = r / scale * Math.cos(lat) * Math.cos(lon); // x
      pos[3 * i + 1] = r / scale * Math.cos(lat) * Math.sin(lon); // y
      pos[3 * i + 2] = r / scale * Math.sin(lat); // x

      // rotation NED (graphics)
      rotationNED.set(
        -Math.sin(lat) * Math.cos(lon), -Math.sin(lat) * Math.sin(lon), Math.cos(lat), 0,
        -Math.sin(lon), Math.cos(lon), 0, 0,
        -Math.cos(lat) * Math.cos(lon), -Math.cos(lat) * Math.sin(lon), -Math.sin(lat), 0,
        0, 0, 0, 1,
      );
      rotationNED.transpose();
      quaternionNED.setFromRotationMatrix(rotationNED);
      quaternionNED.toArray(rotN, i * 4);

      // rotation BODY (graphics)
      rotB[4 * i] = e1;
      rotB[4 * i + 1] = e2;
      rotB[4 * i + 2] = e3;
      rotB[4 * i + 3] = e0;
    }
    // check for NaNs and exlude them from the arrays
    // filter the results and set them to the object arr property
    arr.t = time.filter(el => !Number.isNaN(el));
    arr.alt = alt.filter(el => !Number.isNaN(el));
    arr.pos = pos.filter(el => !Number.isNaN(el));
    arr.rotN = rotN.filter(el => !Number.isNaN(el));
    arr.rotB = rotB.filter(el => !Number.isNaN(el));
    arr.xVel = xVel.filter(el => !Number.isNaN(el));
    arr.yVel = yVel.filter(el => !Number.isNaN(el));
    arr.zVel = zVel.filter(el => !Number.isNaN(el));
    arr.vel = vel.filter(el => !Number.isNaN(el));
    arr.qInf = qInf.filter(el => !Number.isNaN(el));
    arr.tempNose = tempNose.filter(el => !Number.isNaN(el));
    arr.AoA = AoA.filter(el => !Number.isNaN(el));
  },

  /**
  * Creates numeric.spline objects that can be used for further interpolation.
  * Saves calculation.
  *
  * @param {object} arr - arrays with results
  * @param {object} spline - empty object for the results
  *
  */
  createSplines(arr, spline) {
    spline.tVStempNose = numeric.spline(arr.t, arr.tempNose);
    spline.tVSqInf = numeric.spline(arr.t, arr.qInf);
    spline.tVSvel = numeric.spline(arr.t, arr.vel);
    spline.tVSxVel = numeric.spline(arr.t, arr.xVel);
    spline.tVSyVel = numeric.spline(arr.t, arr.yVel);
    spline.tVSzVel = numeric.spline(arr.t, arr.zVel);
    spline.tVSalt = numeric.spline(arr.t, arr.alt);
    spline.tVSAoA = numeric.spline(arr.t, arr.AoA);
  },
};


// ////////////////////////
// ///////////////////////////////// VIEW /////////////////////////////////////////
// ////////////////////////

// ///////////////////////
// LIVE PLOTS (MODULE 4)//
// ///////////////////////
/**
* This module uses Chart.js to define live plots on the website
*/
const livePlots = (function () {
  // creates an array of colours that are used for consecutive datasets (new plots)
  const colours = ['#d67ad5', '#8481e6', '#81cee6', '#9ae681', '#e6e481', '#e6a381'];
  // minimum altitude line
  const dataset1 = {
    label: 'Minimum altitude line',
    data: [{ x: 0, y: 8 }, { x: 8, y: 8 }],
    borderColor: '#db2525',
    borderWidth: 2,
    pointRadius: 0,
    hitRadius: 0,
    pointHoverRadius: 0,
    showLine: true,
  };
  // altitude vs velocity
  const dataset2 = {
    colourIndex: 0, // changes before another plot is plotted on top
    label: 'Altitude vs Velocity',
    data: [], // initialize to being empty
    pointRadius: 0.9,
    get pointBackgroundColor() { return colours[this.colourIndex]; },
    pointHoverRadius: 2,
    pointHoverBackgroundColor: '#FFF',
  };
  // angle of attack vs time
  const dataset3 = {
    colourIndex: 0, // changes before another plot is plotted on top
    label: 'AoA vs time',
    data: [], // initialize to being empty
    pointRadius: 0.8,
    get pointBackgroundColor() { return colours[this.colourIndex]; },
    pointHoverRadius: 1.5,
    pointHoverBackgroundColor: '#FFF',
  };
  // puts altitude vs velocity and minimum altitude line together
  function data1() {
    return { datasets: [dataset1, dataset2] };
  }
  // initializes data for the AoA vs time plot
  function data2() {
    return { datasets: [dataset3] };
  }
  const options1 = {
    // define initial options
    layout: { padding: { left: 10, right: 25, top: 25, bottom: 5 } }, // padding of the graph
    showLines: false, // no lines for the plot
    animation: { duration: 0 }, // general animation time
    hover: {
      mode: 'nearest', // sets which elements appear in the tooltip
      intersect: false, // if true, the hover mode only applies when the mouse position intersects an item on the chart
      animationDuration: 0, // duration in milliseconds it takes to animate hover style changes
    },
    responsive: false, // the chart doesnt resize
    // define event functions
    onClick(evt, point) {
      // test if results.arr object is empty
      if (Object.keys(results.arr).length !== 0 && results.arr.constructor === Object) {
        const index = point[0]._index;
        const altitude = dataset2.data[index].y * 1000; // [m] it's monotonous unlike velocity array
        let time = everpolate.linear(altitude, results.arr.alt, results.arr.t)[0];
        time = Math.round(time * 100) / 100; // round to 2 decimal places
        G.setTime(time); // sets the time of the simulation
      }
    },
    scales: {
      // xAxes configuration
      xAxes: [{
        scaleLabel: {
          display: true,
          labelString: 'Velocity [km/s]',
        },
        gridLines: { color: 'rgba(255, 255, 255, 0.4)' },
        ticks: {
          fontSize: 9,
          stepSize: 1,
          beginAtZero: true, // include 0 on x axis
        },
      }],
      // yAxes configuration
      yAxes: [{
        scaleLabel: {
          display: true,
          labelString: 'Altitude [km]',
        },
        gridLines: { color: 'rgba(255, 255, 255, 0.4)' },
        ticks: {
          fontSize: 9,
          beginAtZero: true, // include 0 on y axis
        },
      }],
    },
    // legend configuration
    legend: { display: false },
    // tooltip configuration
    tooltips: {
      mode: 'nearest',
      intersect: false,
      borderWidth: 0.8,
      borderColor: '#FFF',
    },
  };
  const options2 = {
    // define initial options
    layout: { padding: { left: 10, right: 25, top: 25, bottom: 5 } }, // padding of the graph
    showLines: false,
    animation: { duration: 0 }, // general animation time
    hover: {
      mode: 'nearest',
      intersect: false,
      animationDuration: 0, // duration of animations when hovering an item
    },
    responsive: false, // the chart doesnt resize
    scales: {
      // xAxes configuration
      xAxes: [{
        scaleLabel: {
          display: true,
          labelString: 'Time [s]',
        },
        gridLines: { color: 'rgba(255, 255, 255, 0.4)' },
        ticks: {
          fontSize: 9,
          stepSize: 10,
          beginAtZero: true, // include 0 on x axis
        },
      }],
      // yAxes configuration
      yAxes: [{
        scaleLabel: {
          display: true,
          labelString: 'Angle of Attack [deg]',
        },
        gridLines: { color: 'rgba(255, 255, 255, 0.4)' },
        ticks: {
          fontSize: 9,
          beginAtZero: true, // include 0 on y axis
        },
      }],
    },
    // legend configuration
    legend: { display: false },
    // tooltip configuration
    tooltips: {
      mode: 'nearest',
      intersect: false,
      borderWidth: 0.8,
      borderColor: '#FFF',
    },
  };
  let myPlot1 = null;
  let myPlot2 = null;
  /**
  * Initializes charts in canvas elements
  */
  function init() {
    // gets DOM elements to place the graph
    const ctx1 = document.querySelector('#plot1');
    const ctx2 = document.querySelector('#plot2');

    // changes default font options
    Chart.defaults.global.defaultFontColor = '#FFF';
    Chart.defaults.global.defaultFontFamily = 'Verdana';
    Chart.defaults.global.defaultFontSize = '13';

    // creates 2 plots with predefined data and options
    myPlot1 = new Chart(ctx1, { type: 'scatter', data: data1(), options: options1 }); // altitude vs velocity
    myPlot2 = new Chart(ctx2, { type: 'scatter', data: data2(), options: options2 }); // angle of attack vs time
  }
  /**
  * Updates data on the plots in the render loop
  *
  * @param {number} time - real simulation time
  * @param {object} spline - splines from result object
  */
  function update(time, { spline }) {
    // PLOT 1 (altitude vs velocity)

    // interpolates values in the results arrays to find the current values of altitude and velocity
    // and rounds the values
    const altitude = Math.round(spline.tVSalt.at(time)) / 1000; // [km]
    const velocity = Math.round(spline.tVSvel.at(time)) / 1000; // [km/s]

    const last1 = myPlot1.data.datasets.length - 1; // index of the last item in the datasets array
    if (myPlot1.data.datasets[last1].data.length !== 0) { // if the last dataset is not empty
      const previous1 = myPlot1.data.datasets[last1].data[myPlot1.data.datasets[last1].data.length - 1]; // defines the last data point
      const distance1 = Math.sqrt(400 * (velocity - previous1.x) * (velocity - previous1.x) + (altitude - previous1.y) * (altitude - previous1.y)); // defines the distance between the current and last data point

      // based on distance either push to the last dataset or create a new dataset
      if (distance1 > 1 && distance1 < 50) {
        myPlot1.data.datasets[last1].data.push({ x: velocity, y: altitude }); // push data to the last dataset
        myPlot1.update(0); // triggers an update of the chart. This can be safely called after updating the data object. This will update all scales, legends, and then re-render the chart
      }
      if (distance1 > 50) { // this big a distance signifies a new plot
        const colourIndex = myPlot1.data.datasets[last1].colourIndex + 1; // changes the colour of the next plot
        // initializes new dataset
        myPlot1.data.datasets.push({
          colourIndex, // changes before another plot is plotted on top
          label: 'Altitude vs Velocity',
          data: [],
          pointRadius: 0.9,
          get pointBackgroundColor() { return colours[this.colourIndex]; },
          pointHoverRadius: 2,
          pointHoverBackgroundColor: '#FFF',
        });
        myPlot1.data.datasets[myPlot1.data.datasets.length - 1].data.push({ x: velocity, y: altitude }); // push data to the last dataset
        myPlot1.update(0); // triggers an update of the chart
      }
    } else { // if the last dataset is empty
      myPlot1.data.datasets[last1].data.push({ x: velocity, y: altitude }); // push data to the last dataset
      myPlot1.update(0); // triggers an update of the chart
    }

    // PLOT 2 (angle of attack vs time)

    // interpolates values in the results arrays to find the current values of angle of attack and time
    // and rounds the values
    const t = Math.round(time * 100) / 100; // seconds, round to 2 decimal places
    const AoA = Math.round(spline.tVSAoA.at(time) * 180 / Math.PI * 1000) / 1000; // deg, round to 3 decimal places

    const last2 = myPlot2.data.datasets.length - 1; // index of the last item in the datasets array
    if (myPlot2.data.datasets[last2].data.length !== 0) { // if the last dataset is not empty
      const previous2 = myPlot2.data.datasets[last2].data[myPlot2.data.datasets[last2].data.length - 1]; // defines the last data point
      const distance2 = Math.sqrt((t - previous2.x) * (t - previous2.x) + (AoA - previous2.y) * (AoA - previous2.y)); // defines distance between the current and the last data point

      // based on distance either push to the last dataset or create a new dataset
      if (distance2 > 0.2 && distance2 < 50) {
        myPlot2.data.datasets[last2].data.push({ x: t, y: AoA });
        myPlot2.update(0);
      }

      if (distance2 > 50) { // this big a distance signifies a new plot
        const colourIndex = myPlot2.data.datasets[last2].colourIndex + 1; // changes the colour of the next plot
        // initializes new dataset
        myPlot2.data.datasets.push({
          colourIndex, // changes before another plot is plotted on top
          label: 'AoA vs time',
          data: [{ x: 0, y: trajectoryOptions.cond.alpha * 180 / Math.PI }],
          pointRadius: 0.8,
          get pointBackgroundColor() { return colours[this.colourIndex]; },
          pointHoverRadius: 1.5,
          pointHoverBackgroundColor: '#FFF',
        });
        myPlot2.data.datasets[myPlot2.data.datasets.length - 1].data.push({ x: t, y: AoA }); // push data to the last dataset
        myPlot2.update(0); // triggers an update of the chart
      }
    } else { // if the last dataset is empty
      myPlot2.data.datasets[last2].data.push({ x: t, y: AoA }); // push data to the last dataset
      myPlot2.update(0); // triggers an update of the chart
    }
  }

  // returns public API
  return {
    init,
    update,
  };
})();


// /////////////////////////
// LIVE OUTPUTS (MODULE 5)//
// /////////////////////////
/**
* This module defines live outputs on the website.
*/
const liveOutputs = (function () {
  // grabs DOM elements used to display live outputs
  const DOMel = {
    alt: document.querySelector('#altitude'),
    vel: document.querySelector('#velocity'),
    qInf: document.querySelector('#dynamicPress'),
    tempNose: document.querySelector('#tempNose'),
  };
  // initialize live outputs based on the results from trajectory module
  function init() {
    // defines initial parameters
    const altitude = Math.round(trajectoryOptions.cond.h);
    const velocity = Math.round(trajectoryOptions.cond.vInf);
    const [, , , , , , , , , dynamicPress, mach, , ,] = trajectoryOptions.calcAtmoProperties(altitude, velocity, lander); // eslint-disable-line comma-spacing

    // initializes those values in the live outputs
    DOMel.alt.innerHTML = `${altitude} m`;
    DOMel.vel.innerHTML = `${velocity} m/s`;
    DOMel.qInf.innerHTML = '0 Pa';
    DOMel.tempNose.innerHTML = '0 K';

    // check if the current values of those parameters are within limits to safely deploy a parachute
    checkLimits(altitude, mach, dynamicPress);
  }
  /**
  * Updates altitude, velocity and dynamic pressure outputs in the render loop and checks for acceptable limits of those values.
  *
  * @param {number} time - real simulation time
  * @param {object} spline - splines from result object
  */
  function update(time, { spline }) {
    // gets current values of parameters
    const altitude = Math.round(spline.tVSalt.at(time));
    const velocity = Math.round(spline.tVSvel.at(time));
    const dynamicPress = Math.round(spline.tVSqInf.at(time) * 100) / 100;
    const tempNose = Math.round(spline.tVStempNose.at(time) * 100) / 100;
    const mach = trajectoryOptions.calcAtmoProperties(altitude, velocity, lander)[10];

    // updates values in the DOM elements
    DOMel.alt.innerHTML = `${altitude} m`;
    DOMel.vel.innerHTML = `${velocity} m/s`;
    DOMel.qInf.innerHTML = `${dynamicPress} Pa`;
    DOMel.tempNose.innerHTML = `${tempNose} K`;

    // check if the current values of those parameters are within limits to safely deploy a parachute
    checkLimits(altitude, mach, dynamicPress);
  }
  /**
  * This function checks limits of the current live outputs and colours the font accordingly and displays a message
  * @param {number} altitude - current altitude
  * @param {number} mach - current Mach number
  * @param {number} dynamicPressure - current dynamic pressure
  */
  function checkLimits(altitude, mach, dynamicPressure) {
    // defines limits
    const altitudeLimit = 8000; // [m]
    const machMax = 2.2; // [-]
    const machMin = 1.1; // [-]
    const dynPressMax = 850; // [Pa]
    const dynPressMin = 239; // [Pa]

    // defines conditions for the values
    const altCond = altitude > altitudeLimit;
    const machCond = (mach < machMax && mach > machMin);
    const dynPressCond = (dynamicPressure < dynPressMax && dynamicPressure > dynPressMin);

    // checks for the limits of the numbers and colours them accordingly
    altCond ? DOMel.alt.style.color = 'green' : DOMel.alt.style.color = 'red'; // altitude limits
    machCond ? DOMel.vel.style.color = 'green' : DOMel.vel.style.color = 'red'; // mach number limits
    dynPressCond ? DOMel.qInf.style.color = 'green' : DOMel.qInf.style.color = 'red'; // dynamic pressure limits

    // grabs success and fail message containers
    const successMessage = document.querySelector('#successMessage');
    const failureMessage = document.querySelector('#failureMessage');

    // checks for meeting the right conditions to deploy a parachute
    if (altCond && machCond && dynPressCond) {
      successMessage.style.visibility = 'visible'; // display the success message
    }

    // checks for conditions after crossing which one cannot deploy a parachute
    if (!altCond) {
      failureMessage.style.visibility = 'visible'; // display the failure message
    }
  }

  // returns public API
  return {
    init,
    update,
  };
})();

// ///////////////////////////////////
// GRAPHICS AND ANIMATION (MODULE 6)//
// ///////////////////////////////////

const G = (function () {
  // initializes main variables needed to create 3D graphics
  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.01, 1000);
  camera.position.set(-11.822327, -3.95953, 32.913791);
  camera.up = new THREE.Vector3(0, 0, 1);

  // let temp = new THREE.Vector3();
  const renderer = new THREE.WebGLRenderer({ antialias: true });

  // configures renderer and adds it to the DOM
  const container = document.querySelector('#webgl');
  renderer.setSize(container.clientWidth, container.clientHeight);
  renderer.setPixelRatio(window.devicePixelRatio);
  container.appendChild(renderer.domElement);
  // turn on the physically correct lighting model
  renderer.physicallyCorrectLights = true;
  // configures full screen options
  THREEx.FullScreen.bindKey({ charCode: 'm'.charCodeAt(0) });

  // creates an object for keeping track of time
  const clock = new THREE.Clock();

  // creates and docks an object for monitoring performance
  const stats = new Stats();
  document.body.appendChild(stats.dom);

  // creates a loading manager (handles and keeps track of loaded and pending data)
  const manager = new THREE.LoadingManager();
  const progressBarElem = document.querySelector('#loading .progressbar'); // grabs a progress bar element from DOM
  manager.onProgress = (url, itemsLoaded, itemsTotal) => { // defines behaviour for the loading page
    const progress = itemsLoaded / itemsTotal;
    progressBarElem.style.transform = `scaleX(${progress})`;
  };


  // adds your own window resize functionality
  function onWindowResize() {
    // sets the aspect ratio to match the new browser window aspect ratio
    camera.aspect = container.clientWidth / container.clientHeight;
    // updates the camera's frustum
    camera.updateProjectionMatrix();
    // updates the size of the renderer AND the canvas
    renderer.setSize(container.clientWidth, container.clientHeight);
  }
  window.addEventListener('resize', onWindowResize);

  // //////////////////////////////////////
  //         MESHES IN THE SCENE
  // //////////////////////////////////////

  /**
   * Creates a starry sky in the scene
   */
  const SkyObj = (function (loadingManager) {
    const amount = 3000; // defines the number of stars
    // creates a stars group
    const stars = new THREE.Group();
    stars.name = 'stars';
    // creates a star for each loop
    let size; let geometry; let star;
    // defines material for a star
    const material = new THREE.MeshBasicMaterial({ color: 0xffffff });
    for (let i = 0; i < amount; i++) {
      // makes a star
      size = Math.random() * 0.18 + 0.08; // randomize the size of each star, between 0.08 and 0.26
      geometry = new THREE.SphereBufferGeometry(size, 8, 8);
      star = new THREE.Mesh(geometry, material);

      const radius = 100; // the radius of the sphere on which the stars will appear
      const z = Math.random() * 2 * radius - radius; // randomizes the z-coordinate
      const lam = Math.random() * 2 * Math.PI; // randomizes the longitude
      // calculates x and y coordinates from the known variables
      const x = radius * Math.cos(lam) * Math.sqrt(1 - (z / radius) * (z / radius));
      const y = radius * Math.sin(lam) * Math.sqrt(1 - (z / radius) * (z / radius));
      // assigns positions to a star
      star.position.z = z;
      star.position.x = x;
      star.position.y = y;

      // adds a star to the stars group
      stars.add(star);
      // Three.js Cleanup
      geometry.dispose();
    }
    // Three.js Cleanup
    material.dispose();
    // adds all the stars to the scene
    scene.add(stars);
  })(manager);

  /**
   * Creates a Mars object in the scene
   */
  const MarsObj = (function (loadingManager) {
    // grabs scale and radius of the Mars
    const { scale } = lander.mesh;
    const radius = Mars.r / scale;
    // console.log(radius);
    // creates Mars coordinate system
    const marsFrame = createAxes('PLANET', radius + 10); // + 10);

    // MATERIAL //
    // defines the scale for the bump map (perceived depth of the surface)
    const bumpScale = 0.5;
    // creates a loader for textures
    const loader = new THREE.TextureLoader(manager);
    // initializes texture variables and an array for Mars materials
    let textureMap; let textureBumpMap; let displacementMap;
    const materialsArray = [];
    // for each sphere element
    for (let i = 0; i < 8; i++) {
      // creates textures for the colour map and bump map
      textureMap = loader.load(`images/color/mars${i}.png`);
      textureBumpMap = loader.load(`images/bump/mars${i}.png`);
      displacementMap = textureBumpMap;
      // sets encoding for both maps
      textureMap.encoding = THREE.LinearEncoding; // default
      textureBumpMap.encoding = THREE.LinearEncoding; // default
      // reduces texture blurring by setting the anistropic filtering level
      // only needs to be done for the colour texture
      textureMap.anisotropy = renderer.capabilities.getMaxAnisotropy();
      // assigns textures to the material
      materialsArray[i] = new THREE.MeshPhongMaterial({
        map: textureMap,
        bumpMap: textureBumpMap,
        bumpScale,
        displacementMap,
        displacementScale: 0.8,
        displacementBias: -0.28, // set this so that the mars radius is correct (look at coordinates)
      });
    }

    // GEOMETRY //
    // creates all 8 sphere parts
    const mars = new THREE.Group(); // creates a container for all sphere parts
    let material; let geometry; let mesh; // define variables used for each part

    for (let i = 0; i < 8; i++) {
      geometry = new THREE.SphereBufferGeometry(radius, 60, 60, Math.PI / 2 * (i % 4), Math.PI / 2, Math.PI / 2 * Math.floor(i / 4), Math.PI / 2);
      material = materialsArray[i];
      // makes the surface of the Mars matt
      material.specular = new THREE.Color('black');
      material.shininess = 0;
      // creates a mesh
      mesh = new THREE.Mesh(geometry, material);
      mesh.frustumCulled = false; // all mars parts get rendered every frame even if it isn't visible, this is to avoid scene lags
      mars.add(mesh);
      // THREE.js Cleanup
      materialsArray[i].dispose();
      material.dispose();
      geometry.dispose();
    }

    // TRANSFORMATIONS //
    // adds mars group to the planet frame
    marsFrame.add(mars);
    // to apply following matrix directly, also the object (mars + coordinate system) is static - no need to recalculate matrix
    marsFrame.traverse((obj) => { obj.matrixAutoUpdate = false; });
    // transformation of Mars mesh to match planet frame coordinate system
    mars.matrix.makeRotationFromEuler(new THREE.Euler(Math.PI / 2, Math.PI, 0));
    // adds mars to the scene
    scene.add(marsFrame);

    // returns marsFrame objects
    return marsFrame;
  })(manager);

  // const MarsAtmo = (function () {
  //   // create atmosphere
  //   // grabs scale and radius of the Mars
  //   const { scale } = lander.mesh;
  //   const radius = Mars.r / scale;
  //
  //   const geometry = new THREE.SphereGeometry(1.01 * radius, 64, 64);
  //   const material = THREEx.createAtmosphereMaterial();
  //   material.uniforms.glowColor.value.set(0xfda600);
  //   material.uniforms.coeficient.value = 0.8;
  //   material.uniforms.power.value = 2.2;
  //   // material.side = THREE.DoubleSide;
  //   const mesh = new THREE.Mesh(geometry, material);
  //   mesh.scale.multiplyScalar(1.01);
  //
  //   scene.add(mesh);
  // })();

  /**
   * Creates a lander object in the scene and defines methods on the lander
   */
  const LanderObj = (function (loadingManager) {
    // defines scaling and the curve for the lander
    let { scale, landerMagnification, rLath, xLath } = lander.mesh;
    // defines position and orientation in space
    let { lat, lon, r, e0, e1, e2, e3 } = trajectoryOptions.cond;

    // defines THREE.js constants
    const materials = [new THREE.MeshBasicMaterial({
      color: 0x505050,
      wireframe: false,
      transparent: true,
      side: THREE.DoubleSide,
      // map: loader.load('images/carbon-fiber-texture.jpg'),
    }), new THREE.MeshBasicMaterial({
      color: 0xffffff,
      wireframe: true,
      wireframeLinecap: 'round', // default
      wireframeLinejoin: 'round', // default
      wireframeLinewidth: 0.5,
    })];

    /**
    * Updates variables that influence the lander object and were
    * changed by the GUI
    * @param {string} - specifies what type of update is needed (mesh, position, orientation)
    */
    function update(type) {
      // updates parameters for all types of updates
      ({ scale, landerMagnification, rLath, xLath } = lander.mesh);
      ({ lat, lon, r, e0, e1, e2, e3 } = trajectoryOptions.cond);

      // gets references to the BODY and NED frames
      const BODYframe = scene.getObjectByName('BODY');
      const NEDframe = scene.getObjectByName('NED');

      // main code, executes based on the type of the update
      switch (type) {
        case 'mesh': {
          // gets reference to the current mesh of the lander
          const oldvehicle = scene.getObjectByName('MESH');
          // removes old lander
          BODYframe.remove(oldvehicle);
          // creates new vehicle geometry, material and mesh
          const points = [];
          for (let i = 0; i < xLath.length; i++) {
            points.push(new THREE.Vector2(rLath[i] / scale * landerMagnification, xLath[i] / scale * landerMagnification));
          }
          // const loader = new THREE.TextureLoader();
          const geometry = new THREE.LatheBufferGeometry(points, rLath.length - 1); // a shape generated by spinning a line
          const vehicle = SceneUtils.createMultiMaterialObject(geometry, materials);
          vehicle.name = 'MESH';
          // makes parent child relationships between the objects
          BODYframe.add(vehicle);
          // sets initial rotation of the mesh wrt the BODY frame
          vehicle.rotation.z = Math.PI / 2;
          // Three.js cleanup
          geometry.dispose();
          break;
        }
        case 'position': {
          // calculates new planet frame coordinates of the lander
          const x = r / scale * Math.cos(lat) * Math.cos(lon);
          const y = r / scale * Math.cos(lat) * Math.sin(lon);
          const z = r / scale * Math.sin(lat);
          // sets the new position
          NEDframe.position.set(x, y, z);
          break;
        }
        case 'orientation': {
          // updates BODY orientation
          BODYframe.quaternion.set(e1, e2, e3, e0); // (x, y, z, w)
          // updates NED orientation
          const rotation = new THREE.Matrix4();
          // rotation matrix (https://en.wikipedia.org/wiki/North_east_down)
          rotation.set(
            -Math.sin(lat) * Math.cos(lon), -Math.sin(lat) * Math.sin(lon), Math.cos(lat), 0,
            -Math.sin(lon), Math.cos(lon), 0, 0,
            -Math.cos(lat) * Math.cos(lon), -Math.cos(lat) * Math.sin(lon), -Math.sin(lat), 0,
            0, 0, 0, 1,
          );
          rotation.transpose();
          NEDframe.quaternion.setFromRotationMatrix(rotation);
          break;
        }
        default:
      }
    }
    /**
    * Creates a lander mesh in the initial position and orientation.
    * Creates two coordinate frames: NED and BODY frames.
    *
    */
    function init() {
      // creates mesh group for BODY and NED FRAMES with axes helpers
      const BODYframe = createAxes('BODY', 0.15);
      const NEDframe = createAxes('NED', 25);

      // creates vehicle geometry and mesh
      const points = [];
      for (let i = 0; i < xLath.length; i++) {
        points.push(new THREE.Vector2(rLath[i] / scale * landerMagnification, xLath[i] / scale * landerMagnification));
      }
      // const loader = new THREE.TextureLoader();
      const geometry = new THREE.LatheBufferGeometry(points, rLath.length - 1); // a shape generated by spinning a line
      const vehicle = SceneUtils.createMultiMaterialObject(geometry, materials);
      vehicle.name = 'MESH';
      // makes parent-child relationships between the objects
      BODYframe.add(vehicle);
      NEDframe.add(BODYframe);
      scene.add(NEDframe);

      // ROTATIONS //
      // set initial rotation of the mesh wrt the BODY frame
      vehicle.rotation.z = Math.PI / 2;
      // set initial rotation of the BODY frame wrt the NED frame (quaternions)
      BODYframe.quaternion.set(e1, e2, e3, e0); // (x, y, z, w)
      // set initial rotation of the NED frame wrt the PLANET frame (latitude and longitude)
      const rotation = new THREE.Matrix4();
      // rotation matrix (https://en.wikipedia.org/wiki/North_east_down)
      rotation.set(
        -Math.sin(lat) * Math.cos(lon), -Math.sin(lat) * Math.sin(lon), Math.cos(lat), 0,
        -Math.sin(lon), Math.cos(lon), 0, 0,
        -Math.cos(lat) * Math.cos(lon), -Math.cos(lat) * Math.sin(lon), -Math.sin(lat), 0,
        0, 0, 0, 1,
      );
      rotation.transpose();
      NEDframe.quaternion.setFromRotationMatrix(rotation);

      // POSITIONS //
      // set initial position of the NED frame (latitude, longitude and altitude) in PLANET FRAME
      // using properties
      const xPLANET = r / scale * Math.cos(lat) * Math.cos(lon);
      const yPLANET = r / scale * Math.cos(lat) * Math.sin(lon);
      const zPLANET = r / scale * Math.sin(lat);
      NEDframe.position.set(xPLANET, yPLANET, zPLANET);

      /**
       * Creates light in the scene which is added to the Lander object so that it can follow it
       */
      const light = new THREE.DirectionalLight(0xffffff, 4);
      light.position.set(0, 0, -10);
      // const helper = new THREE.DirectionalLightHelper(light, 5);
      // scene.add(helper);
      NEDframe.add(light);

      // THREE.js cleanup
      geometry.dispose();
    }
    return {
      init,
      update,
    };
  })(manager);
  LanderObj.init();


  /**
   * Creates a velocity object in the scene and defines its methods
   */
  const Velocity = (function () {
    /**
    * creates initial velocity vector
    */
    function init() {
      // defines velocity coordinates
      const { u, v, w, vInf } = trajectoryOptions.cond;
      // grabs a NED frame from the scene
      const NEDframe = scene.getObjectByName('NED');

      // defines direction, origin and length of the velocity vector
      const dir = new THREE.Vector3(u, v, w);
      dir.normalize(); // has to be normalized
      const origin = new THREE.Vector3(0, 0, 0);
      const length = vInf / 100000;

      // creates velocity vector and assign a name to it
      const arrowHelper = new THREE.ArrowHelper(dir, origin, length, 0xffffff, 0.006, 0.002);
      arrowHelper.name = 'VELOCITY';

      // adds the vector to NEDframe
      NEDframe.add(arrowHelper);
    }
    /**
    * updates velocity vector based on the GUI input or in the animation loop
    * @param {string} type - "GUI" or "animation"
    * @param {number} time - current time in the render loop
    */
    function update(type, time) {
      // gets reference to the velocity object
      const velVector = scene.getObjectByName('VELOCITY');

      // based on whether we are updating the velocity in the GUI or in the render loop, defines velocity components
      let u; let v; let w; let vInf;
      switch (type) {
        case 'GUI': {
          // updates values of velocity components
          ({ u, v, w, vInf } = trajectoryOptions.cond);
          break;
        }
        case 'animation': {
          // grabs spline results for the update in the render loop
          const { spline } = results;
          // interpolates values in the results arrays to find current velocity values
          u = spline.tVSxVel.at(time);
          v = spline.tVSyVel.at(time);
          w = spline.tVSzVel.at(time);
          vInf = spline.tVSvel.at(time);
          break;
        }
        default:
      }
      // calculates new direction and length of the velocity vector
      const length = vInf / 100000;
      const dir = new THREE.Vector3(u, v, w);
      dir.normalize();
      // sets new direction and length of the vector in the scene
      velVector.setDirection(dir);
      velVector.setLength(length, 0.006, 0.003);
    }


    return {
      init,
      update,
    };
  })();
  Velocity.init();

  /**
   * Configures camera and controls settings dependent on the NED frame.
   */
  const controls = new THREE.TrackballControls(camera, container);
  function config() {
    // grabs a NED frame
    const NEDframe = scene.getObjectByName('NED');
    // const BODYframe = scene.getObjectByName('BODY');

    // positions camera
    // and configures trackball controlls

    controls.update();

    // configures controls
    controls.target = NEDframe.position;
    controls.rotateSpeed = 2;
    controls.zoomSpeed = 2;
    controls.noPan = true;
    controls.minDistance = 0.1;
    controls.maxDistance = 95;
    controls.staticMoving = true;
  }
  config();
  /**
  * Creates an axes helper representing frame of reference.
  * Used in createMars and createLander functions.
  *
  * @param {string} name name property assigned to the group object
  * @param {number} axesLength length of the coordinate system axis
  *
  * @return {object} THREE.js group object representing the frame of reference.
  */
  function createAxes(name, axesLength) {
    // creates a new group for the frame
    const frame = new THREE.Group();
    // assigns a name to the group
    frame.name = name;
    // creates axes and adds it to the group
    frame.add(new THREE.AxesHelper(axesLength));

    return frame;
  }

  // //////////////////////////////////////
  //        ANIMATION SYSTEM
  // //////////////////////////////////////
  // defines animation variables
  let mixer; let action; let isPlay = false;
  /**
   * Defines keyframes, clip, mixer, and clipAction objects for the animation system.
   *
   * @param {object} results.arr - arrays with the results
   */
  function defineAnim({ arr }) {
    // grabs NED frame object
    const NEDframe = scene.getObjectByName('NED');

    // creates a keyframes track for NED ROTATION
    const quaternionNKF = new THREE.QuaternionKeyframeTrack('.quaternion', arr.t, arr.rotN);

    // creates a keyframes track for NED POSITION
    const positionKF = new THREE.VectorKeyframeTrack('.position', arr.t, arr.pos);

    // creates a keyframes track for BODY ROTATION
    const quaternionBKF = new THREE.QuaternionKeyframeTrack('.children[1].quaternion', arr.t, arr.rotB); // BODYframe is children[1] of NEDframe

    // creates an Animation Clip from the tracks
    const tracks = [quaternionNKF, positionKF, quaternionBKF];
    const duration = -1; // use -1 to automatically calculate the length from the array of tracks
    const clip = new THREE.AnimationClip('entry', duration, tracks);
    clip.optimize(); // optimizes each track by removing equivalent sequential keys

    // sets up the Animation Mixer for the lander
    mixer = new THREE.AnimationMixer(NEDframe);

    // creates an Animation Action and configures it
    action = mixer.clipAction(clip);
    action.loop = THREE.LoopOnce; // makes the animation play once

    // adds event listener for when the action finishes
    mixer.addEventListener('finished', (e) => {
      action.stop(); // this method also calls reset()
      isPlay = false; // the animation is no longer playing
      Velocity.update('GUI'); // updates the velocity vector
      liveOutputs.init(); // updates the live outputs
    });
  }

  // play animation
  function playAnim() {
    // does not respond when there are no results
    if (Object.keys(results.arr).length === 0) {
      const calcnotpressed = document.querySelector('#calcnotpressed');
      calcnotpressed.style.visibility = 'visible';
      return;
    }

    if (!isPlay) {
      action.play();
      action.timeScale = 1;
      isPlay = true;
      action.paused = false;
    } else { // if the animation is already playing only change the speed
      action.timeScale = 1;
    }
  }

  // plays animation 5 times faster
  function playAnimX5() {
    // does not respond when there are no results
    if (Object.keys(results.arr).length === 0) {
      const calcnotpressed = document.querySelector('#calcnotpressed');
      calcnotpressed.style.visibility = 'visible';
      return;
    }

    if (!isPlay) {
      action.play();
      action.timeScale = 5;
      isPlay = true;
      action.paused = false;
    } else { // if the animation is already playing only change the speed
      action.timeScale = 5;
    }
  }

  // pauses running animation
  function pauseAnim() {
    // does not respond when there are no results
    if (Object.keys(results.arr).length === 0) {
      const calcnotpressed = document.querySelector('#calcnotpressed');
      calcnotpressed.style.visibility = 'visible';
      return;
    }

    action.paused = true;
    isPlay = false;
  }

  // stops running animation
  function stopAnim() {
    // does not respond when there are no results
    if (Object.keys(results.arr).length === 0) {
      const calcnotpressed = document.querySelector('#calcnotpressed');
      calcnotpressed.style.visibility = 'visible';
      return;
    }

    action.stop(); // this method also calls reset()
    isPlay = false;
  }

  // sets clip Action time
  function setTime(time) { action.time = time; }

  /**
   * Calls requestAnimationFrame in a loop and renders the scene.
   */
  function animate() {
    updateView();
    renderer.render(scene, camera); // renders the scene
    window.requestAnimationFrame(animate); // calls the animate function in a loop
  }

  /**
   * Renders the scene and updates changing variables in the scene.
   */
  function updateView() {
    const delta = clock.getDelta(); // get time passed from last callback

    if (mixer) {
      const currentTime = action.time; // get the current real simulation time
      // if animation is not playing stop updating live outputs
      if (isPlay) {
        Velocity.update('animation', currentTime); // updates velocity vector
        livePlots.update(currentTime, results); // updates live plots
        liveOutputs.update(currentTime, results); // updates live outputs
      }
      mixer.update(delta);
    }

    // define the behaviour of the camera when it is about to cross the Mars surface
    const currentPosition = camera.position; // the current camera position (three js vector)
    const distance = currentPosition.distanceTo(new THREE.Vector3(0, 0, 0)); // calculate the distance between the centre of Mars and the camera
    if (distance < 34) { // if the distance is smaller than Mars radius then change the position of the camera
      const { x, y, z } = currentPosition; // current camera coordinates
      const delr = 34 - distance; // desired change of distance
      const sign = (x > 0) ? 1 : -1;
      // based on your calculations
      const delx = sign * 2 * delr / Math.sqrt(1 + (y / x) * (y / x) + (z / x) * (z / x));
      const dely = y * delx / x;
      const delz = z * delx / x;
      // set camera position
      camera.position.set(x + delx, y + dely, z + delz);
      camera.up = new THREE.Vector3(0, 0, 1);
    }


    controls.update(); // updates camera controls
    stats.update(); // updates the performance monitor
  }

  // return public API
  return {
    camera,
    scene,
    manager,

    MarsObj,
    LanderObj,
    Velocity,

    animate,
    setTime,
    defineAnim,
    playAnim,
    playAnimX5,
    pauseAnim,
    stopAnim,
  };
})();

// /////////////////////////////////////
// GRAPHICAL USER INTERFACE (MODULE 7)//
// /////////////////////////////////////

// expand this module to be a mediator (an object literal)- a single point of communication between all other modules
// use dat.gui to achieve it

/**
 * Creates GUI with initial parameters.
 *
 * @param {object} lander - contains geometry and lander properties
 * @param {object} trajectoryOptions - contains atmosphere data
 */
function setupGUI() {
  // initializes GUI and sets its width
  const gui = new dat.GUI({ width: 400 }); // specify width of gui

  const { cond } = trajectoryOptions; // contains initial conditions, it is the same object, its just a pointer

  // add save menu, remember the initial state of objects before we implement any changes
  gui.remember(lander);
  gui.remember(cond);

  // MASS of the heatshield
  gui.add(lander, 'mass', 10, 50000).name('Mass [kg]');
  // VELOCITY of the heatshield
  const vel = gui.add(cond, 'vInf', 1000, 10000).name('Velocity [m/s]');
  vel.onChange(() => { G.Velocity.update('GUI'); liveOutputs.init(); });

  // GEOMETRY of the heatshield
  let folder = gui.addFolder('GEOMETRY');
  // number of mesh divisions
  const nop = folder.add(lander.mesh, 'nop', 10, 50).step(1).name('Mesh Divisions');
  nop.onChange(() => { lander.defineMesh(); G.LanderObj.update('mesh'); });
  // diameter of the heatshield
  const dia = folder.add(lander, 'diameter', 1, 30).name('Diameter [m]');
  dia.onChange(() => { lander.defineMesh(); G.LanderObj.update('mesh'); });
  // nose radius
  const nRad = folder.add(lander, 'rNose', 0.1, 10).name('Nose Radius [m]');
  nRad.onChange(() => { lander.defineMesh(); G.LanderObj.update('mesh'); });
  // angle
  const theta = folder.add(lander, 'sphereConeAngleDeg', 10, 80).name('Cone Angle [deg]');
  theta.onChange(() => { lander.defineMesh(); G.LanderObj.update('mesh'); });

  // POSITION
  folder = gui.addFolder('POSITION');
  // altitude
  const alt = folder.add(cond, 'h', 100000, 200000).step(5000).name('Altitude [m]');
  alt.onChange(() => { G.LanderObj.update('position'); liveOutputs.init(); });
  // latitude
  const lat = folder.add(cond, 'latDeg', -90, 90).step(5).name('Latitude [deg]');
  lat.onChange(() => { G.LanderObj.update('position'); G.LanderObj.update('orientation'); });
  // longitude
  const lon = folder.add(cond, 'lonDeg', 0, 360).step(5).name('Longitude [deg]');
  lon.onChange(() => { G.LanderObj.update('position'); G.LanderObj.update('orientation'); });

  // ANGLES
  folder = gui.addFolder('ANGLES');
  // flight path angle
  const FPA = folder.add(cond, 'FPADeg', 0, 90).name('FPA [deg]');
  FPA.onChange(() => { G.Velocity.update('GUI'); G.LanderObj.update('orientation'); });
  // azimuth
  const azi = folder.add(cond, 'aziDeg', -180, 180).name('Azimuth [deg]');
  azi.onChange(() => { G.Velocity.update('GUI'); G.LanderObj.update('orientation'); });
  // angle of attack
  const AoA = folder.add(cond, 'alphaDeg', -80, 80).name('Angle of Attack [deg]');
  AoA.onChange(() => G.LanderObj.update('orientation'));

  // ATMOSPHERES
  gui.add(trajectoryOptions, 'current', { Pathfinder: 'pathfinder', Curiosity: 'curiosity', Opportunity: 'opportunity', Phoenix: 'phoenix', Schiaparelli: 'schiaparelli', Spirit: 'spirit' }).name('Atmosphere');

  // TRAJECTORY function - triggers calculation of the trajectory
  const run = gui.add(results, 'run').name('CALCULATE');

  // create folder for animation
  folder = gui.addFolder('ANIMATION');
  // play
  folder.add(G, 'playAnim').name('Play');
  // play x5
  folder.add(G, 'playAnimX5').name('Play Faster');
  // stop
  const stop = folder.add(G, 'stopAnim').name('Stop');
  stop.onFinishChange(() => G.Velocity.update('GUI'));
  // pause
  folder.add(G, 'pauseAnim').name('Pause');
}

// /////////////////////////
// LOGIC
// /////////////////////////

// /////////////////////////
// ////////////////////////END * OF * THE * SCRIPT///////////////////////////////////
// /////////////////////////
