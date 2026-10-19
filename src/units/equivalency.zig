//! Physical equivalencies: conversions between units of *different* dimensions that
//! are related by a physical law or convention (frequency <-> wavelength, mass <->
//! energy, redshift <-> velocity, and so on). Ordinary same-dimension conversions
//! are handled by `Quantity.to`; these are the cross-dimension cases that astropy
//! exposes through its `equivalencies` mechanism.
//!
//! Each equivalency is a namespaced struct exposing a `convert` function. Most use a
//! hub-and-spoke design: the input is first converted to a single canonical unit
//! (the "hub"), then from the hub to the requested output. This keeps the number of
//! branches linear in the number of supported units rather than quadratic. The
//! `Equivalency` enum at the bottom ties them together for use from
//! `Quantity.toWithEquivalency`.
//!
//! Equivalencies that need extra information (a rest frequency, a beam area, a
//! zero-point flux, ...) take a comptime-checked `args` tuple; those that don't
//! ignore it. Argument validation lives in each struct's `checkArgs`, which emits a
//! `@compileError` on misuse.
//!

const std = @import("std");
const testing = std.testing;
const dim_mod = @import("dim.zig");
const Dim = dim_mod.Dim;
const unit_mod = @import("unit.zig");
const Unit = unit_mod.Unit;
const UNITLESS = unit_mod.UNITLESS;
const quantity_mod = @import("quantity.zig");
const Quantity = quantity_mod.Quantity;
const utils = @import("utils.zig");
const si = @import("si.zig");
const photo = @import("photometric.zig");
const misc = @import("misc.zig");
const cst = @import("../constants.zig");
const imperial = @import("imperial.zig");
const us = @import("units_storage.zig");
const zatest = @import("../test.zig");

fn unknownConversion(comptime name: []const u8) []const u8 {
    return "Unknown " ++ name ++ " equivalency conversion.";
}

/// Spectral equivalency: converts between the various ways of locating a point in a
/// spectrum -- frequency, wavelength, energy, spectroscopic wavenumber (m-1), and
/// angular wavenumber (rad/m). Hub unit is Hz.
///
/// Relations (with c the speed of light, h the Planck constant):
///   nu = c / lambda,  E = h nu,  nu = c * wavenumber_spec,  nu = c * k / (2 pi)
pub const spectral_equivalency = struct {
    const name = "spectral";
    const Hz: Unit = si.Hz;
    const dim_freq: Dim = Hz.dim;
    const dim_length: Dim = si.m.dim;
    const dim_energy: Dim = si.J.dim;
    const dim_wavenumber_spec: Dim = si.m.pow(-1).dim;
    const dim_wavenumber_ang: Dim = si.rad.div(si.m).dim;
    const two_pi: Quantity(f64, si.rad) = .init(2 * std.math.pi);
    const h = cst.h.quantity;
    const c = cst.c.quantity;

    pub fn convert(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit) Quantity(T, U_out) {
        const as_Hz: Quantity(T, Hz) = toHz(T, U_in, input_quantity);
        return fromHz(T, &as_Hz, U_out);
    }

    pub fn convertInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, out: *Quantity(T, U_out)) void {
        var temp_Hz: Quantity(T, si.Hz) = .init(out.value);
        toHzInto(T, U_in, input_quantity, &temp_Hz);
        fromHzInto(T, &temp_Hz, U_out, out);
    }

    fn toHz(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in)) Quantity(T, si.Hz) {
        if (comptime U_in.dim.eql(dim_freq)) {
            return input_quantity.to(Hz);
        } else if (comptime U_in.dim.eql(dim_length)) {
            return input_quantity.pow(-1).mulScalar(c).to(Hz);
        } else if (comptime U_in.dim.eql(dim_energy)) {
            return input_quantity.divScalar(h).to(Hz);
        } else if (comptime U_in.dim.eql(dim_wavenumber_spec)) {
            return input_quantity.mulScalar(c).to(Hz);
        } else if (comptime U_in.dim.eql(dim_wavenumber_ang)) {
            return input_quantity.mulScalar(c).divScalar(two_pi).to(Hz);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn toHzInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), out: *Quantity(T, Hz)) void {
        if (comptime U_in.dim.eql(dim_freq)) {
            input_quantity.toInto(Hz, out);
        } else if (comptime U_in.dim.eql(dim_length)) {
            var temp_pow: Quantity(T, U_in.pow(-1)) = .init(out.value);
            input_quantity.powInto(-1, &temp_pow);
            var temp_mul_c: Quantity(T, @TypeOf(temp_pow).unit.mul(@TypeOf(c).unit)) = .init(temp_pow.value);
            temp_pow.mulScalarInto(c, &temp_mul_c);
            temp_mul_c.toInto(Hz, out);
        } else if (comptime U_in.dim.eql(dim_energy)) {
            var temp_div_h: Quantity(T, U_in.div(@TypeOf(h).unit)) = .init(out.value);
            input_quantity.divScalarInto(h, &temp_div_h);
            temp_div_h.toInto(Hz, out);
        } else if (comptime U_in.dim.eql(dim_wavenumber_spec)) {
            var temp_mul_c: Quantity(T, U_in.mul(@TypeOf(c).unit)) = .init(out.value);
            input_quantity.mulScalarInto(c, &temp_mul_c);
            temp_mul_c.toInto(Hz, out);
        } else if (comptime U_in.dim.eql(dim_wavenumber_ang)) {
            const c_div_two_pi = c.div(two_pi);
            var temp_mul: Quantity(T, U_in.mul(@TypeOf(c_div_two_pi))) = .init(out.value);
            input_quantity.mulScalarInto(c_div_two_pi, &temp_mul);
            temp_mul.toInto(Hz, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromHz(comptime T: type, input_quantity: *const Quantity(T, Hz), comptime U_out: Unit) Quantity(T, U_out) {
        if (comptime U_out.dim.eql(dim_freq)) {
            return input_quantity.to(U_out);
        } else if (comptime U_out.dim.eql(dim_length)) {
            return input_quantity.pow(-1).mulScalar(c).to(U_out);
        } else if (comptime U_out.dim.eql(dim_energy)) {
            return input_quantity.mulScalar(h).to(U_out);
        } else if (comptime U_out.dim.eql(dim_wavenumber_spec)) {
            return input_quantity.divScalar(c).to(U_out);
        } else if (comptime U_out.dim.eql(dim_wavenumber_ang)) {
            return input_quantity.mulScalar(two_pi).divScalar(c).to(U_out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromHzInto(comptime T: type, input_quantity: *const Quantity(T, Hz), comptime U_out: Unit, out: *Quantity(T, U_out)) void {
        if (comptime U_out.dim.eql(dim_freq)) {
            input_quantity.toInto(U_out, out);
        } else if (comptime U_out.dim.eql(dim_length)) {
            var temp_pow: Quantity(T, Hz.pow(-1)) = .init(out.value);
            input_quantity.powInto(-1, &temp_pow);
            var temp_mul_c: Quantity(T, Hz.pow(-1).mul(@TypeOf(c).unit)) = .init(out.value);
            temp_pow.mulScalarInto(c, &temp_mul_c);
            temp_mul_c.toInto(U_out, out);
        } else if (comptime U_out.dim.eql(dim_energy)) {
            var temp_mul_h: Quantity(T, Hz.mul(@TypeOf(h).unit)) = .init(out.value);
            input_quantity.mulScalarInto(h, &temp_mul_h);
            temp_mul_h.toInto(U_out, out);
        } else if (comptime U_out.dim.eql(dim_wavenumber_spec)) {
            var temp_div_c: Quantity(T, Hz.div(@TypeOf(c))) = .init(out.value);
            input_quantity.divScalarInto(c, &temp_div_c);
            temp_div_c.toInto(U_out, out);
        } else if (comptime U_out.dim.eql(dim_wavenumber_ang)) {
            const two_pi_div_c = two_pi.div(c);
            var temp_mul: Quantity(T, Hz.mul(@TypeOf(two_pi_div_c).unit)) = .init(out.value);
            input_quantity.mulScalarInto(two_pi_div_c, &temp_mul);
            temp_mul.toInto(U_out, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }
};

/// Mass-energy equivalency: converts between mass and energy via E = m c^2.
/// Hub unit is the joule.
// TODO: Support generalise equivalency as Astropy does such as kg m-2 <-> J m-2 (??).
pub const mass_energy_equivalency = struct {
    const name = "mass-energy";
    const J: Unit = si.J;
    const dim_energy: Dim = si.J.dim;
    const dim_mass: Dim = si.kg.dim;
    const c = cst.c.quantity;

    pub fn convert(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit) Quantity(T, U_out) {
        const as_J: Quantity(T, J) = toJ(T, U_in, input_quantity);
        return fromJ(T, &as_J, U_out);
    }

    pub fn convertInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, out: *Quantity(T, U_out)) void {
        var temp: Quantity(T, J) = .init(out.value);
        toJInto(T, U_in, input_quantity, &temp);
        fromJInto(T, &temp, U_out, out);
    }

    fn toJ(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in)) Quantity(T, J) {
        if (comptime U_in.dim.eql(dim_energy)) {
            return input_quantity.to(J);
        } else if (comptime U_in.dim.eql(dim_mass)) {
            return input_quantity.mulScalar(c.pow(2)).to(J);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn toJInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), out: *Quantity(T, J)) void {
        if (comptime U_in.dim.eql(dim_energy)) {
            input_quantity.toInto(J, out);
        } else if (comptime U_in.dim.eql(dim_mass)) {
            const c2 = c.pow(2);
            var temp_mul: Quantity(T, U_in.mul(@TypeOf(c2).unit)) = .init(out.value);
            input_quantity.mulScalarInto(c2, &temp_mul);
            temp_mul.toInto(J, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromJ(comptime T: type, input_quantity: *const Quantity(T, J), comptime U_out: Unit) Quantity(T, U_out) {
        if (comptime U_out.dim.eql(dim_energy)) {
            return input_quantity.to(U_out);
        } else if (comptime U_out.dim.eql(dim_mass)) {
            return input_quantity.divScalar(c.pow(2)).to(U_out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromJInto(comptime T: type, input_quantity: *const Quantity(T, J), comptime U_out: Unit, out: *Quantity(T, U_out)) void {
        if (comptime U_out.dim.eql(dim_energy)) {
            input_quantity.toInto(U_out, out);
        } else if (comptime U_out.dim.eql(dim_mass)) {
            const c2 = c.pow(2);
            var temp: Quantity(T, J.div(@TypeOf(c2).unit)) = .init(out.value);
            input_quantity.divScalarInto(c2, &temp);
            temp.toInto(U_out, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }
};

/// Temperature-energy equivalency: converts between thermal energy and temperature
/// via E = k_B T, where k_B is the Boltzmann constant. Hub unit is the joule.
/// Temperature inputs are first normalized to kelvin so offset scales (e.g. Celsius)
/// are handled correctly before the physics is applied.
pub const temperature_energy_equivalency = struct {
    const name = "temperature-energy";
    const J: Unit = si.J;
    const K: Unit = si.K;
    const dim_energy: Dim = si.J.dim;
    const dim_temp: Dim = si.K.dim;
    const k_B = cst.k_B.quantity;

    pub fn convert(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit) Quantity(T, U_out) {
        const as_J: Quantity(T, J) = toJ(T, U_in, input_quantity);
        return fromJ(T, &as_J, U_out);
    }

    pub fn convertInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, out: *Quantity(T, U_out)) void {
        var temp: Quantity(T, J) = .init(out.value);
        toJInto(T, U_in, input_quantity, &temp);
        fromJInto(T, &temp, U_out, out);
    }

    fn toJ(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in)) Quantity(T, J) {
        if (comptime U_in.dim.eql(dim_energy)) {
            return input_quantity.to(J);
        } else if (comptime U_in.dim.eql(dim_temp)) {
            return input_quantity.to(si.K).mulScalar(k_B).to(J);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn toJInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), out: *Quantity(T, J)) void {
        if (comptime U_in.dim.eql(dim_energy)) {
            input_quantity.toInto(J, out);
        } else if (comptime U_in.dim.eql(dim_temp)) {
            var temp_to_K: Quantity(T, K) = .init(out.value);
            input_quantity.toInto(K, &temp_to_K);
            var temp_mul: Quantity(T, K.mul(@TypeOf(k_B).unit)) = .init(temp_to_K.value);
            temp_to_K.mulScalarInto(k_B, &temp_mul);
            temp_mul.toInto(J, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromJ(comptime T: type, input_quantity: *const Quantity(T, J), comptime U_out: Unit) Quantity(T, U_out) {
        if (comptime U_out.dim.eql(dim_energy)) {
            return input_quantity.to(U_out);
        } else if (comptime U_out.dim.eql(dim_temp)) {
            return input_quantity.divScalar(k_B).to(U_out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromJInto(comptime T: type, input_quantity: *const Quantity(T, J), comptime U_out: Unit, out: *Quantity(T, U_out)) void {
        if (comptime U_out.dim.eql(dim_energy)) {
            input_quantity.toInto(U_out, out);
        } else if (comptime U_out.dim.eql(dim_temp)) {
            var temp_div: Quantity(T, J.div(@TypeOf(k_B).unit)) = .init(out.value);
            input_quantity.divScalarInto(k_B, &temp_div);
            temp_div.toInto(U_out, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }
};

// In astropy, this equivalency is required because offset-based units (°C, °F)
// are not handled by regular .to(). In this implementation, offset handling is
// built into the unit system, so this equivalency simply delegates to .to().
/// Temperature equivalency: provided for API parity with astropy, where it is needed
/// because offset scales (Celsius, Fahrenheit) are otherwise incompatible. Here unit
/// offsets are built into `Quantity.to`, so this simply delegates to it.
pub const temperature_equivalency = struct {
    pub fn convert(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit) Quantity(T, U_out) {
        return input_quantity.to(U_out);
    }

    pub fn convertInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, out: *Quantity(T, U_out)) void {
        input_quantity.toInto(U_out, out);
    }
};

/// Molar-mass / atomic-mass-unit equivalency: converts between a molar mass (g/mol)
/// and an atomic mass (u), related by Avogadro's number N_A (1 u corresponds to
/// 1 g/mol). Hub unit is the atomic mass unit.
pub const molar_mass_amu_equivalency = struct {
    const name = "mola-mass-amu";
    const molar_mass: Unit = si.g.div(si.mol);
    const molar_mass_dim: Dim = molar_mass.dim;
    const amu: Unit = misc.u;
    const amu_dim: Dim = amu.dim;
    const N_A = cst.N_A.quantity;

    pub fn convert(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit) Quantity(T, U_out) {
        const as_amu: Quantity(T, amu) = toAmu(T, U_in, input_quantity);
        return fromAmu(T, &as_amu, U_out);
    }

    pub fn convertInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, out: *Quantity(T, U_out)) void {
        var temp_amu: Quantity(T, amu) = .init(out.value);
        toAmuInto(T, U_in, input_quantity, &temp_amu);
        fromAmuInto(T, &temp_amu, U_out, out);
    }

    fn toAmu(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in)) Quantity(T, amu) {
        if (comptime U_in.dim.eql(amu_dim)) {
            return input_quantity.to(amu);
        } else if (comptime U_in.dim.eql(molar_mass_dim)) {
            return input_quantity.divScalar(N_A).to(amu);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn toAmuInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), out: *Quantity(T, amu)) void {
        if (comptime U_in.dim.eql(amu_dim)) {
            input_quantity.toInto(amu, out);
        } else if (comptime U_in.dim.eql(molar_mass_dim)) {
            var temp_div: Quantity(T, U_in.div(@TypeOf(N_A).unit)) = .init(out.value);
            input_quantity.divScalarInto(N_A, &temp_div);
            temp_div.toInto(amu, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromAmu(comptime T: type, input_quantity: *const Quantity(T, amu), comptime U_out: Unit) Quantity(T, U_out) {
        if (comptime U_out.dim.eql(amu_dim)) {
            return input_quantity.to(U_out);
        } else if (comptime U_out.dim.eql(molar_mass_dim)) {
            return input_quantity.mulScalar(N_A).to(U_out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromAmuInto(comptime T: type, input_quantity: *const Quantity(T, amu), comptime U_out: Unit, out: *Quantity(T, U_out)) void {
        if (comptime U_out.dim.eql(amu_dim)) {
            input_quantity.toInto(U_out, out);
        } else if (comptime U_out.dim.eql(molar_mass_dim)) {
            var temp_mul: Quantity(T, amu.mul(@TypeOf(N_A).unit)) = .init(out.value);
            input_quantity.mulScalarInto(N_A, &temp_mul);
            temp_mul.toInto(U_out, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }
};

/// Doppler-redshift equivalency: converts between dimensionless redshift z and radial
/// velocity using the relativistic relation. Hub unit is km/s.
///   z -> v:  v = c * ((z+1)^2 - 1) / ((z+1)^2 + 1)
///   v -> z:  z = sqrt((1 + beta)/(1 - beta)) - 1,  with beta = v/c
pub const doppler_redshift_equivalency = struct {
    const name = "doppler-redshift";
    const rv_unit: Unit = us.km.div(si.s);
    const rv_dim: Dim = rv_unit.dim;
    const c_kms = cst.c.quantity.to(rv_unit);

    pub fn convert(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit) Quantity(T, U_out) {
        const as_rv: Quantity(T, rv_unit) = toRv(T, U_in, input_quantity);
        return fromRv(T, &as_rv, U_out);
    }

    pub fn convertInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, out: *Quantity(T, U_out)) void {
        var temp: Quantity(T, rv_unit) = .init(out.value);
        toRvInto(T, U_in, input_quantity, &temp);
        fromRvInto(T, &temp, U_out, out);
    }

    fn toRv(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in)) Quantity(T, rv_unit) {
        if (comptime U_in.dim.eql(rv_dim)) {
            return input_quantity.to(rv_unit);
        } else if (comptime U_in.dim.eql(UNITLESS.dim)) {
            const zponesq = input_quantity.addScalarValue(1).pow(2);
            return zponesq.subScalarValue(1).div(zponesq.addScalarValue(1)).mulScalar(c_kms).to(rv_unit);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn toRvInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), out: *Quantity(T, rv_unit)) void {
        if (comptime U_in.dim.eql(rv_dim)) {
            input_quantity.toInto(rv_unit, out);
        } else if (comptime U_in.dim.eql(UNITLESS.dim)) {
            var buf: Quantity(T, UNITLESS) = .init(out.value);
            input_quantity.toInto(UNITLESS, &buf); // buf = z
            buf.addScalarValueInPlace(1); // buf = z + 1
            buf.powInto(2, &buf); // buf = (z + 1)^2 = x
            buf.addScalarValueInPlace(1); // buf = x + 1
            buf.powInto(-1, &buf); // buf = 1 / (x + 1)
            buf.mulScalarValueInPlace(-2); // buf = -2 / (x + 1)
            buf.addScalarValueInPlace(1); // buf = 1 - 2 / (x + 1) = (x - 1)/(x + 1)

            var vel: Quantity(T, UNITLESS.mul(rv_unit)) = .init(out.value);
            buf.mulScalarInto(c_kms, &vel); // vel = c . (x - 1)/(x + 1) = v
            vel.toInto(rv_unit, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromRv(comptime T: type, input_quantity: *const Quantity(T, rv_unit), comptime U_out: Unit) Quantity(T, U_out) {
        if (comptime U_out.dim.eql(rv_dim)) {
            return input_quantity.to(U_out);
        } else if (comptime U_out.dim.eql(UNITLESS.dim)) {
            const beta = input_quantity.divScalar(c_kms).to(UNITLESS);
            return beta.addScalarValue(1).divScalar(beta.neg().addScalarValue(1)).sqrt().subScalarValue(1).to(U_out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromRvInto(comptime T: type, input_quantity: *const Quantity(T, rv_unit), comptime U_out: Unit, out: *Quantity(T, U_out)) void {
        if (comptime U_out.dim.eql(rv_dim)) {
            input_quantity.toInto(U_out, out);
        } else if (comptime U_out.dim.eql(UNITLESS.dim)) {
            var beta_raw: Quantity(T, rv_unit.div(rv_unit)) = .init(out.value);
            input_quantity.divScalarInto(c_kms, &beta_raw); // v/c, unit rv/rv
            var buf: Quantity(T, UNITLESS) = .init(out.value);
            beta_raw.toInto(UNITLESS, &buf); // beta
            buf.subScalarValueInPlace(1); // beta - 1
            buf.powInto(-1, &buf); // 1 / (beta - 1)
            buf.mulScalarValueInPlace(-2); // -2 / (beta - 1) = 2 / (1 - beta)
            buf.subScalarValueInPlace(1); // 2 / (1 - beta) - 1 = (1 + beta)/(1 - beta)
            buf.sqrtInto(&buf); // sqrt(...)
            buf.subScalarValueInPlace(1); // z
            buf.toInto(U_out, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }
};

/// Dimensionless-angles equivalency: lets angular quantities (radians, degrees, ...)
/// be treated as dimensionless and vice versa, normalizing through radians. Useful
/// because angle is dimensionless in SI but tracked as a dimension here. Hub is the
/// dimensionless value (with the angle stripped to / restored from radians).
pub const dimensionless_angles_equivalency = struct {
    const name = "dimensionless-angles";

    pub fn convert(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit) Quantity(T, U_out) {
        const as_unitless: Quantity(T, UNITLESS) = toUnitless(T, U_in, input_quantity);
        return fromUnitless(T, &as_unitless, U_out);
    }

    pub fn convertInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, out: *Quantity(T, U_out)) void {
        var temp: Quantity(T, UNITLESS) = .init(out.value);
        toUnitlessInto(T, U_in, input_quantity, &temp);
        fromUnitlessInto(T, &temp, U_out, out);
    }

    fn toUnitless(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in)) Quantity(T, UNITLESS) {
        const u_in_dim: Dim = U_in.dim;
        if (comptime U_in.eqlExact(UNITLESS)) {
            return input_quantity.*;
        } else if (comptime u_in_dim.l.eqlScalar(0) and
            u_in_dim.m.eqlScalar(0) and
            u_in_dim.t.eqlScalar(0) and
            u_in_dim.i.eqlScalar(0) and
            u_in_dim.th.eqlScalar(0) and
            u_in_dim.n.eqlScalar(0) and
            u_in_dim.j.eqlScalar(0) and
            !u_in_dim.a.eqlScalar(0))
        {
            return .init(input_quantity.to(si.rad).value);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn toUnitlessInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), out: *Quantity(T, UNITLESS)) void {
        const u_in_dim: Dim = U_in.dim;
        if (comptime U_in.eqlExact(UNITLESS)) {
            @memcpy(out.value, input_quantity.value);
        } else if (comptime u_in_dim.l.eqlScalar(0) and
            u_in_dim.m.eqlScalar(0) and
            u_in_dim.t.eqlScalar(0) and
            u_in_dim.i.eqlScalar(0) and
            u_in_dim.th.eqlScalar(0) and
            u_in_dim.n.eqlScalar(0) and
            u_in_dim.j.eqlScalar(0) and
            !u_in_dim.a.eqlScalar(0))
        {
            var temp: Quantity(T, si.rad) = .init(out.value);
            input_quantity.toInto(si.rad, &temp);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromUnitless(comptime T: type, input_quantity: *const Quantity(T, UNITLESS), comptime U_out: Unit) Quantity(T, U_out) {
        const u_out_dim: Dim = U_out.dim;
        if (comptime U_out.eqlExact(UNITLESS)) {
            return input_quantity.*;
        } else if (comptime u_out_dim.l.eqlScalar(0) and
            u_out_dim.m.eqlScalar(0) and
            u_out_dim.t.eqlScalar(0) and
            u_out_dim.i.eqlScalar(0) and
            u_out_dim.th.eqlScalar(0) and
            u_out_dim.n.eqlScalar(0) and
            u_out_dim.j.eqlScalar(0) and
            !u_out_dim.a.eqlScalar(0))
        {
            return Quantity(T, si.rad).init(input_quantity.value).to(U_out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromUnitlessInto(comptime T: type, input_quantity: *const Quantity(T, UNITLESS), comptime U_out: Unit, out: *Quantity(T, U_out)) void {
        const u_out_dim: Dim = U_out.dim;
        if (comptime U_out.eqlExact(UNITLESS)) {
            // No-op: the sole caller (convertInto) has already written the result
            // into out's buffer, which input_quantity aliases. UNITLESS<-UNITLESS
            // is the identity, so there is nothing to copy.
        } else if (comptime u_out_dim.l.eqlScalar(0) and
            u_out_dim.m.eqlScalar(0) and
            u_out_dim.t.eqlScalar(0) and
            u_out_dim.i.eqlScalar(0) and
            u_out_dim.th.eqlScalar(0) and
            u_out_dim.n.eqlScalar(0) and
            u_out_dim.j.eqlScalar(0) and
            !u_out_dim.a.eqlScalar(0))
        {
            var temp: Quantity(T, si.rad) = .init(input_quantity.value);
            temp.toInto(U_out, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }
};

// NOTE: At the moment, contrary to Astropy, this conversion gives negative angles and negative
// distances, instead of NaN.
/// Parallax equivalency: converts between a parallax angle and a distance via the
/// reciprocal relation p[arcsec] = 1 / d[pc]. Hub unit is the arcsecond.
pub const parallax_equivalency = struct {
    const name = "parallax";
    const as: Unit = si.arcsec;
    const as_dim: Dim = as.dim;
    const pc: Unit = us.pc;
    const pc_dim: Dim = pc.dim;

    pub fn convert(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit) Quantity(T, U_out) {
        const as_arcsec: Quantity(T, si.arcsec) = toArcsec(T, U_in, input_quantity);
        return fromArcsec(T, &as_arcsec, U_out);
    }

    pub fn convertInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, out: *Quantity(T, U_out)) void {
        var temp: Quantity(T, si.arcsec) = .init(out.value);
        toArcsecInto(T, U_in, input_quantity, &temp);
        fromArcsecInto(T, &temp, U_out, out);
    }

    fn toArcsec(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in)) Quantity(T, as) {
        if (comptime U_in.dim.eql(as_dim)) {
            return input_quantity.to(as);
        } else if (comptime U_in.dim.eql(pc_dim)) {
            return .init(input_quantity.to(pc).pow(-1).value);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn toArcsecInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), out: *Quantity(T, as)) void {
        if (comptime U_in.dim.eql(as_dim)) {
            input_quantity.toInto(as, out);
        } else if (comptime U_in.dim.eql(pc_dim)) {
            var temp: Quantity(T, pc) = .init(out.value);
            input_quantity.toInto(pc, &temp);
            var temp_inv: Quantity(T, pc.pow(-1)) = .init(temp.value);
            temp.powInto(-1, &temp_inv);
            // out shares this buffer; leaving 1/d here is p[arcsec] by definition of parsec.
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromArcsec(comptime T: type, input_quantity: *const Quantity(T, as), comptime U_out: Unit) Quantity(T, U_out) {
        if (comptime U_out.dim.eql(as_dim)) {
            return input_quantity.to(U_out);
        } else if (comptime U_out.dim.eql(pc_dim)) {
            return Quantity(T, pc).init(input_quantity.pow(-1).value).to(U_out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromArcsecInto(comptime T: type, input_quantity: *const Quantity(T, as), comptime U_out: Unit, out: *Quantity(T, U_out)) void {
        if (comptime U_out.dim.eql(as_dim)) {
            input_quantity.toInto(U_out, out);
        } else if (comptime U_out.dim.eql(pc_dim)) {
            var temp_inv: Quantity(T, as.pow(-1)) = .init(out.value);
            input_quantity.powInto(-1, &temp_inv); // out.value = 1/p (numerically d[pc])
            var temp_pc: Quantity(T, pc) = .init(out.value);
            temp_pc.toInto(U_out, out); // reinterpret as pc, then pc -> U_out
        } else {
            @compileError(unknownConversion(name));
        }
    }
};

/// Zero-point-flux equivalency: converts between a linear flux and a "maggy" (a
/// dimensionless flux ratio relative to a reference), via f = F / F0. Requires one
/// argument: the zero-point flux F0 (e.g. 3631 Jy for the AB system). Only two units
/// are involved (flux and dimensionless maggy), so this is a direct conversion rather
/// than hub-and-spoke.
pub const zero_point_flux_equivalency = struct {
    const name = "zero-point-flux";

    fn checkArgs(args: anytype) void {
        const ArgsType = @TypeOf(args);
        const args_fields = @typeInfo(ArgsType).@"struct".fields;

        if (comptime args_fields.len != 1) {
            @compileError("Zero-point-flux equivalency expect 1 argument: the zero-point flux.");
        }
        const ArgType = args_fields[0].type;
        utils.assertIsScalarQuantity(ArgType);

        if (comptime ArgType.unit.dim.eql(UNITLESS.dim)) {
            @compileError("Zero-point-flux equivalency argument must be a flux Quantity, not dimensionless.");
        }
    }

    pub fn convert(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, args: anytype) Quantity(T, U_out) {
        checkArgs(args);
        const flux0 = args.@"0";
        const Flux0Type = @TypeOf(args.@"0");
        const flux0_dim: Dim = comptime Flux0Type.unit.dim;
        const maggy_dim: Dim = UNITLESS.dim;

        if (comptime U_in.dim.eql(maggy_dim) and U_out.dim.eql(flux0_dim)) {
            return input_quantity.mulScalar(flux0).to(U_out);
        } else if (comptime U_in.dim.eql(flux0_dim) and U_out.dim.eql(maggy_dim)) {
            return input_quantity.divScalar(flux0).to(U_out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    pub fn convertInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, out: *Quantity(T, U_out), args: anytype) void {
        checkArgs(args);
        const flux0 = args.@"0";
        const Flux0Type = @TypeOf(args.@"0");
        const flux0_dim: Dim = comptime Flux0Type.unit.dim;
        const maggy_dim: Dim = UNITLESS.dim;

        if (comptime U_in.dim.eql(maggy_dim) and U_out.dim.eql(flux0_dim)) {
            var temp: Quantity(T, U_in.mul(Flux0Type.unit)) = .init(out.value);
            input_quantity.mulScalarInto(flux0, &temp);
            temp.toInto(U_out, out);
        } else if (comptime U_in.dim.eql(flux0_dim) and U_out.dim.eql(maggy_dim)) {
            var temp: Quantity(T, U_in.div(Flux0Type.unit)) = .init(out.value);
            input_quantity.divScalarInto(flux0, &temp);
            temp.toInto(U_out, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }
};

/// Relativistic Doppler equivalency: converts a spectral coordinate (frequency,
/// wavelength, energy, wavenumber) to/from radial velocity using the exact
/// relativistic formula, given a rest spectral location. Requires one argument: the
/// rest value as any spectral-compatible Quantity (normalized to Hz internally).
/// Hub unit is km/s; spectral inputs/outputs route through `spectral_equivalency`.
///   f -> v:  v = c * (f0^2 - f^2) / (f0^2 + f^2)
///   v -> f:  f = f0 * sqrt(1 - beta^2) / (1 + beta),  beta = v/c
pub const doppler_relativistic_equivalency = struct {
    const name = "doppler-relativistic";
    const wavelength_dim: Dim = si.m.dim;
    const energy_dim: Dim = si.J.dim;
    const frequency_dim: Dim = si.Hz.dim;
    const wavenumber_spec_dim: Dim = si.m.pow(-1).dim;
    const wavenumber_ang_dim: Dim = si.rad.div(si.m).dim;
    const rv_unit: Unit = us.km.div(si.s);
    const rv_dim: Dim = rv_unit.dim;
    const c_kms = cst.c.quantity.to(rv_unit);

    fn checkArgs(args: anytype) void {
        const ArgsType = @TypeOf(args);
        const args_fields = @typeInfo(ArgsType).@"struct".fields;

        if (comptime args_fields.len != 1) {
            @compileError("Doppler-relativistic equivalency expect 1 argument.");
        }
        const ArgType = args_fields[0].type;
        utils.assertIsScalarQuantity(ArgType);

        if (comptime !utils.isSpectralUnit(ArgType.unit)) {
            @compileError("Doppler-relativistic equivalency argument must be a spectral compatible Quantity.");
        }
    }

    pub fn convert(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, args: anytype) Quantity(T, U_out) {
        checkArgs(args);
        const rest = args.@"0";
        const RestType = @TypeOf(args.@"0");
        const rest_unit = comptime RestType.unit;
        const rest_hz = spectral_equivalency.convert(RestType.value_type, rest_unit, &rest, si.Hz);
        const value_f64: f64 = @floatCast(rest_hz.value);
        const rest_hz_f64: Quantity(f64, si.Hz) = .init(value_f64);
        const as_rv: Quantity(T, rv_unit) = toRv(T, U_in, input_quantity, &rest_hz_f64);
        return fromRv(T, &as_rv, U_out, &rest_hz_f64);
    }

    pub fn convertInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, out: *Quantity(T, U_out), args: anytype) void {
        checkArgs(args);
        const rest = args.@"0";
        const RestType = @TypeOf(args.@"0");
        const rest_unit = comptime RestType.unit;
        const rest_hz = spectral_equivalency.convert(RestType.value_type, rest_unit, &rest, si.Hz);
        const value_f64: f64 = @floatCast(rest_hz.value);
        const rest_hz_f64: Quantity(f64, si.Hz) = .init(value_f64);
        var temp: Quantity(T, rv_unit) = .init(out.value);
        toRvInto(T, U_in, input_quantity, &rest_hz_f64, &temp);
        fromRvInto(T, &temp, U_out, &rest_hz_f64, out);
    }

    fn toRv(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), rest: *const Quantity(f64, si.Hz)) Quantity(T, rv_unit) {
        if (comptime U_in.dim.eql(rv_dim)) {
            return input_quantity.to(rv_unit);
        } else if (comptime U_in.dim.eql(frequency_dim)) {
            const f = input_quantity.to(si.Hz);
            return fToRv(T, &f, rest);
        } else if (comptime U_in.dim.eql(wavelength_dim) or
            U_in.dim.eql(energy_dim) or
            U_in.dim.eql(wavenumber_spec_dim) or
            U_in.dim.eql(wavenumber_ang_dim))
        {
            const f = spectral_equivalency.convert(T, U_in, input_quantity, si.Hz);
            return fToRv(T, &f, rest);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn toRvInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), rest: *const Quantity(f64, si.Hz), out: *Quantity(T, rv_unit)) void {
        if (comptime U_in.dim.eql(rv_dim)) {
            input_quantity.toInto(rv_unit, out);
        } else if (comptime U_in.dim.eql(frequency_dim)) {
            var temp_f: Quantity(T, si.Hz) = .init(out.value);
            input_quantity.toInto(si.Hz, &temp_f);
            fToRvInto(T, &temp_f, rest, out);
        } else if (comptime U_in.dim.eql(wavelength_dim) or
            U_in.dim.eql(energy_dim) or
            U_in.dim.eql(wavenumber_spec_dim) or
            U_in.dim.eql(wavenumber_ang_dim))
        {
            var temp_f: Quantity(T, si.Hz) = .init(out.value);
            spectral_equivalency.convertInto(T, U_in, input_quantity, si.Hz, &temp_f);
            fToRvInto(T, &temp_f, rest, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromRv(comptime T: type, input_quantity: *const Quantity(T, rv_unit), comptime U_out: Unit, rest: *const Quantity(f64, si.Hz)) Quantity(T, U_out) {
        if (comptime U_out.dim.eql(rv_dim)) {
            return input_quantity.to(U_out);
        } else if (comptime U_out.dim.eql(frequency_dim)) {
            const f: Quantity(T, si.Hz) = rvToF(T, input_quantity, rest);
            return f.to(U_out);
        } else if (comptime U_out.dim.eql(energy_dim) or
            U_out.dim.eql(wavelength_dim) or
            U_out.dim.eql(wavenumber_spec_dim) or
            U_out.dim.eql(wavenumber_ang_dim))
        {
            const f: Quantity(T, si.Hz) = rvToF(T, input_quantity, rest);
            return spectral_equivalency.convert(T, si.Hz, &f, U_out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromRvInto(comptime T: type, input_quantity: *const Quantity(T, rv_unit), comptime U_out: Unit, rest: *const Quantity(f64, si.Hz), out: *Quantity(T, U_out)) void {
        if (comptime U_out.dim.eql(rv_dim)) {
            input_quantity.toInto(U_out, out);
        } else if (comptime U_out.dim.eql(frequency_dim)) {
            var temp_f: Quantity(T, si.Hz) = .init(out.value);
            rvToFInto(T, input_quantity, rest, &temp_f);
            temp_f.toInto(U_out, out);
        } else if (comptime U_out.dim.eql(energy_dim) or
            U_out.dim.eql(wavelength_dim) or
            U_out.dim.eql(wavenumber_spec_dim) or
            U_out.dim.eql(wavenumber_ang_dim))
        {
            var temp_f: Quantity(T, si.Hz) = .init(out.value);
            rvToFInto(T, input_quantity, rest, &temp_f);
            spectral_equivalency.convertInto(T, si.Hz, &temp_f, U_out, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fToRv(comptime T: type, freq: *const Quantity(T, si.Hz), rest: *const Quantity(f64, si.Hz)) Quantity(T, rv_unit) {
        const f0sq = rest.pow(2);
        const fsq = freq.pow(2);
        return (fsq.neg().addScalar(f0sq)).div(fsq.addScalar(f0sq)).mulScalar(c_kms).to(rv_unit);
    }

    fn fToRvInto(comptime T: type, freq: *const Quantity(T, si.Hz), rest: *const Quantity(f64, si.Hz), out: *Quantity(T, rv_unit)) void {
        const f0sq = rest.pow(2);
        const two_f0sq = f0sq.mulScalarValue(2);

        // V/c = 2 f_0^2 / (f^2 + f_0^2) - 1
        var fsq: Quantity(T, si.Hz.pow(2)) = .init(out.value);
        freq.powInto(2, &fsq); // f^2
        fsq.addScalarInPlace(f0sq); // f^2 + f_0^2
        var inv: Quantity(T, si.Hz.pow(-2)) = .init(out.value);
        fsq.powInto(-1, &inv); // 1/(f^2 + f_0^2)
        var ratio: Quantity(T, @TypeOf(inv).unit.mul(@TypeOf(two_f0sq).unit)) = .init(out.value);
        inv.mulScalarInto(two_f0sq, &ratio); // 2f_0^2/(f^2 + f_0^2)
        var ratio_unitless: Quantity([]f64, UNITLESS) = .init(out.value);
        ratio.toInto(UNITLESS, &ratio_unitless);
        ratio_unitless.subScalarValueInPlace(1); // (f_0^2 - f^2)/(f_0^2 + f^2)
        var vel: Quantity(T, UNITLESS.mul(rv_unit)) = .init(out.value);
        ratio_unitless.mulScalarInto(c_kms, &vel); // c . ratio = V
        vel.toInto(rv_unit, out);
    }

    fn rvToF(comptime T: type, rv: *const Quantity(T, rv_unit), rest: *const Quantity(f64, si.Hz)) Quantity(T, si.Hz) {
        const beta = rv.divScalar(c_kms).to(UNITLESS);
        return beta.pow(2).neg().addScalarValue(1).sqrt().div(beta.addScalarValue(1)).mulScalar(rest.*);
    }

    fn rvToFInto(comptime T: type, rv: *const Quantity(T, rv_unit), rest: *const Quantity(f64, si.Hz), out: *Quantity(T, si.Hz)) void {
        const f0 = rest.*;

        // f = f_0 . sqrt(2/(1+beta) - 1), beta = V/c
        var beta_raw: Quantity(T, rv_unit.div(rv_unit)) = .init(out.value);
        rv.divScalarInto(c_kms, &beta_raw); // V/c
        var buf: Quantity(T, UNITLESS) = .init(out.value);
        beta_raw.toInto(UNITLESS, &buf); // beta
        buf.addScalarValueInPlace(1); // 1 + beta
        buf.powInto(-1, &buf); // 1/(1 + beta)
        buf.mulScalarValueInPlace(2); // 2/(1 + beta)
        buf.subScalarValueInPlace(1); // 2/(1 + beta) - 1 = (1 - beta)/(1 + beta)
        buf.sqrtInto(&buf); // sqrt(...)
        var freq: Quantity(T, UNITLESS.mul(si.Hz)) = .init(out.value);
        buf.mulScalarInto(f0, &freq); // f_0 . sqrt(...)
        freq.toInto(si.Hz, out);
    }
};

/// Optical Doppler equivalency: converts a spectral coordinate to/from radial
/// velocity using the optical convention (linear in wavelength). Requires the rest
/// value as a spectral-compatible Quantity. Hub unit is km/s.
///   f -> v:  v = c * (f0 - f) / f
///   v -> f:  f = f0 / (1 + v/c)
pub const doppler_optical_equivalency = struct {
    const name = "doppler-optical";
    const wavelength_dim: Dim = si.m.dim;
    const energy_dim: Dim = si.J.dim;
    const frequency_dim: Dim = si.Hz.dim;
    const wavenumber_spec_dim: Dim = si.m.pow(-1).dim;
    const wavenumber_ang_dim: Dim = si.rad.div(si.m).dim;
    const rv_unit: Unit = us.km.div(si.s);
    const rv_dim: Dim = rv_unit.dim;
    const c_kms = cst.c.quantity.to(rv_unit);

    fn checkArgs(args: anytype) void {
        const ArgsType = @TypeOf(args);
        const args_fields = @typeInfo(ArgsType).@"struct".fields;

        if (comptime args_fields.len != 1) {
            @compileError("Doppler-optical equivalency expect 1 argument.");
        }
        const ArgType = args_fields[0].type;
        utils.assertIsScalarQuantity(ArgType);

        if (comptime !utils.isSpectralUnit(ArgType.unit)) {
            @compileError("Doppler-optical equivalency argument must be a spectral compatible Quantity.");
        }
    }

    pub fn convert(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, args: anytype) Quantity(T, U_out) {
        checkArgs(args);
        const rest = args.@"0";
        const RestType = @TypeOf(args.@"0");
        const rest_unit = comptime RestType.unit;
        const rest_hz = spectral_equivalency.convert(RestType.value_type, rest_unit, &rest, si.Hz);
        const rest_value_f64: f64 = @floatCast(rest_hz.value);
        const rest_hz_f64: Quantity(f64, si.Hz) = .init(rest_value_f64);
        const as_rv: Quantity(T, rv_unit) = toRv(T, U_in, input_quantity, &rest_hz_f64);
        return fromRv(T, &as_rv, U_out, &rest_hz_f64);
    }

    pub fn convertInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, out: *Quantity(T, U_out), args: anytype) void {
        checkArgs(args);
        const rest = args.@"0";
        const RestType = @TypeOf(args.@"0");
        const rest_unit = comptime RestType.unit;
        const rest_hz = spectral_equivalency.convert(RestType.value_type, rest_unit, &rest, si.Hz);
        const rest_value_f64: f64 = @floatCast(rest_hz.value);
        const rest_hz_f64: Quantity(f64, si.Hz) = .init(rest_value_f64);
        var temp: Quantity(T, rv_unit) = .init(out.value);
        toRvInto(T, U_in, input_quantity, &rest_hz_f64, &temp);
        fromRvInto(T, &temp, U_out, &rest_hz_f64, out);
    }

    fn toRv(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), rest: *const Quantity(f64, si.Hz)) Quantity(T, rv_unit) {
        if (comptime U_in.dim.eql(rv_dim)) {
            return input_quantity.to(rv_unit);
        } else if (comptime U_in.dim.eql(frequency_dim)) {
            const f = input_quantity.to(si.Hz);
            return fToRv(T, &f, rest);
        } else if (comptime U_in.dim.eql(wavelength_dim) or
            U_in.dim.eql(energy_dim) or
            U_in.dim.eql(wavenumber_spec_dim) or
            U_in.dim.eql(wavenumber_ang_dim))
        {
            const f = spectral_equivalency.convert(T, U_in, input_quantity, si.Hz);
            return fToRv(T, &f, rest);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn toRvInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), rest: *const Quantity(f64, si.Hz), out: *Quantity(T, rv_unit)) void {
        if (comptime U_in.dim.eql(rv_dim)) {
            input_quantity.toInto(rv_unit, out);
        } else if (comptime U_in.dim.eql(frequency_dim)) {
            var temp_f: Quantity(T, si.Hz) = .init(out.value);
            input_quantity.toInto(si.Hz, &temp_f);
            fToRvInto(T, &temp_f, rest, out);
        } else if (comptime U_in.dim.eql(wavelength_dim) or
            U_in.dim.eql(energy_dim) or
            U_in.dim.eql(wavenumber_spec_dim) or
            U_in.dim.eql(wavenumber_ang_dim))
        {
            var temp_f: Quantity(T, si.Hz) = .init(out.value);
            spectral_equivalency.convertInto(T, U_in, input_quantity, si.Hz, &temp_f);
            fToRvInto(T, &temp_f, rest, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromRv(comptime T: type, input_quantity: *const Quantity(T, rv_unit), comptime U_out: Unit, rest: *const Quantity(f64, si.Hz)) Quantity(T, U_out) {
        if (comptime U_out.dim.eql(rv_dim)) {
            return input_quantity.to(U_out);
        } else if (comptime U_out.dim.eql(frequency_dim)) {
            const f = rvToF(T, input_quantity, rest);
            return f.to(U_out);
        } else if (comptime U_out.dim.eql(energy_dim) or
            U_out.dim.eql(wavelength_dim) or
            U_out.dim.eql(wavenumber_spec_dim) or
            U_out.dim.eql(wavenumber_ang_dim))
        {
            const f = rvToF(T, input_quantity, rest);
            return spectral_equivalency.convert(T, si.Hz, &f, U_out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromRvInto(comptime T: type, input_quantity: *const Quantity(T, rv_unit), comptime U_out: Unit, rest: *const Quantity(f64, si.Hz), out: *Quantity(T, U_out)) void {
        if (comptime U_out.dim.eql(rv_dim)) {
            input_quantity.toInto(U_out, out);
        } else if (comptime U_out.dim.eql(frequency_dim)) {
            var temp_f: Quantity(T, si.Hz) = .init(out.value);
            rvToFInto(T, input_quantity, rest, &temp_f);
            temp_f.toInto(U_out, out);
        } else if (comptime U_out.dim.eql(energy_dim) or
            U_out.dim.eql(wavelength_dim) or
            U_out.dim.eql(wavenumber_spec_dim) or
            U_out.dim.eql(wavenumber_ang_dim))
        {
            var temp_f: Quantity(T, si.Hz) = .init(out.value);
            rvToFInto(T, input_quantity, rest, &temp_f);
            spectral_equivalency.convertInto(T, si.Hz, &temp_f, U_out, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fToRv(comptime T: type, freq: *const Quantity(T, si.Hz), rest: *const Quantity(f64, si.Hz)) Quantity(T, rv_unit) {
        // v = c * (f0 - f) / f
        return (freq.neg().addScalar(rest.*)).div(freq.*).mulScalar(c_kms).to(rv_unit);
    }

    fn fToRvInto(comptime T: type, freq: *const Quantity(T, si.Hz), rest: *const Quantity(f64, si.Hz), out: *Quantity(T, rv_unit)) void {
        const f0 = rest.*;

        // v = c * (f0 - f) / f
        var inv: Quantity(T, si.Hz.pow(-1)) = .init(out.value);
        freq.powInto(-1, &inv); // 1/f
        var ratio: Quantity(T, @TypeOf(inv).unit.mul(@TypeOf(f0).unit)) = .init(out.value);
        inv.mulScalarInto(f0, &ratio); // f_0/f
        var ratio_unitless: Quantity([]f64, UNITLESS) = .init(out.value);
        ratio.toInto(UNITLESS, &ratio_unitless);
        ratio_unitless.subScalarValueInPlace(1); // (f_0 - f)/f
        var vel: Quantity(T, UNITLESS.mul(rv_unit)) = .init(out.value);
        ratio_unitless.mulScalarInto(c_kms, &vel); // c . (f_0 - f)/f = V
        vel.toInto(rv_unit, out);
    }

    fn rvToF(comptime T: type, rv: *const Quantity(T, rv_unit), rest: *const Quantity(f64, si.Hz)) Quantity(T, si.Hz) {
        // f = f0 / (1 + v/c)
        const beta = rv.divScalar(c_kms).to(UNITLESS);
        return (beta.addScalarValue(1)).pow(-1).mulScalar(rest.*);
    }

    fn rvToFInto(comptime T: type, rv: *const Quantity(T, rv_unit), rest: *const Quantity(f64, si.Hz), out: *Quantity(T, si.Hz)) void {
        const f0 = rest.*;

        // f = f0 / (1 + v/c)
        var beta_raw: Quantity(T, rv_unit.div(rv_unit)) = .init(out.value);
        rv.divScalarInto(c_kms, &beta_raw); // V/c
        var buf: Quantity(T, UNITLESS) = .init(out.value);
        beta_raw.toInto(UNITLESS, &buf); // beta
        buf.addScalarValueInPlace(1); // 1 + beta
        buf.powInto(-1, &buf); // 1/(1 + beta)
        var freq: Quantity(T, UNITLESS.mul(si.Hz)) = .init(out.value);
        buf.mulScalarInto(f0, &freq); // f0/(1 + beta)
        freq.toInto(si.Hz, out);
    }
};

/// Radio Doppler equivalency: converts a spectral coordinate to/from radial velocity
/// using the radio convention (linear in frequency). Requires the rest value as a
/// spectral-compatible Quantity. Hub unit is km/s.
///   f -> v:  v = c * (f0 - f) / f0
///   v -> f:  f = f0 * (1 - v/c)
pub const doppler_radio_equivalency = struct {
    const name = "doppler-radio";
    const wavelength_dim: Dim = si.m.dim;
    const energy_dim: Dim = si.J.dim;
    const frequency_dim: Dim = si.Hz.dim;
    const wavenumber_spec_dim: Dim = si.m.pow(-1).dim;
    const wavenumber_ang_dim: Dim = si.rad.div(si.m).dim;
    const rv_unit: Unit = us.km.div(si.s);
    const rv_dim: Dim = rv_unit.dim;
    const c_kms = cst.c.quantity.to(rv_unit);

    fn checkArgs(args: anytype) void {
        const ArgsType = @TypeOf(args);
        const args_fields = @typeInfo(ArgsType).@"struct".fields;

        if (comptime args_fields.len != 1) {
            @compileError("Doppler-radio equivalency expect 1 argument.");
        }
        const ArgType = args_fields[0].type;
        utils.assertIsScalarQuantity(ArgType);

        if (comptime !utils.isSpectralUnit(ArgType.unit)) {
            @compileError("Doppler-radio equivalency argument must be a spectral compatible Quantity.");
        }
    }

    pub fn convert(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, args: anytype) Quantity(T, U_out) {
        checkArgs(args);
        const rest = args.@"0";
        const RestType = @TypeOf(args.@"0");
        const rest_unit = comptime RestType.unit;
        const rest_hz = spectral_equivalency.convert(RestType.value_type, rest_unit, &rest, si.Hz);
        const rest_value_f64: f64 = @floatCast(rest_hz.value);
        const rest_hz_f64: Quantity(f64, si.Hz) = .init(rest_value_f64);
        const as_rv: Quantity(T, rv_unit) = toRv(T, U_in, input_quantity, &rest_hz_f64);
        return fromRv(T, &as_rv, U_out, &rest_hz_f64);
    }

    pub fn convertInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, out: *Quantity(T, U_out), args: anytype) void {
        checkArgs(args);
        const rest = args.@"0";
        const RestType = @TypeOf(args.@"0");
        const rest_unit = comptime RestType.unit;
        const rest_hz = spectral_equivalency.convert(RestType.value_type, rest_unit, &rest, si.Hz);
        const rest_value_f64: f64 = @floatCast(rest_hz.value);
        const rest_hz_f64: Quantity(f64, si.Hz) = .init(rest_value_f64);
        var temp: Quantity(T, rv_unit) = .init(out.value);
        toRvInto(T, U_in, input_quantity, &rest_hz_f64, &temp);
        fromRvInto(T, &temp, U_out, &rest_hz_f64, out);
    }

    fn toRv(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), rest: *const Quantity(f64, si.Hz)) Quantity(T, rv_unit) {
        if (comptime U_in.dim.eql(rv_dim)) {
            return input_quantity.to(rv_unit);
        } else if (comptime U_in.dim.eql(frequency_dim)) {
            const f = input_quantity.to(si.Hz);
            return fToRv(T, &f, rest);
        } else if (comptime U_in.dim.eql(wavelength_dim) or
            U_in.dim.eql(energy_dim) or
            U_in.dim.eql(wavenumber_spec_dim) or
            U_in.dim.eql(wavenumber_ang_dim))
        {
            const f = spectral_equivalency.convert(T, U_in, input_quantity, si.Hz);
            return fToRv(T, &f, rest);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn toRvInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), rest: *const Quantity(f64, si.Hz), out: *Quantity(T, rv_unit)) void {
        if (comptime U_in.dim.eql(rv_dim)) {
            input_quantity.toInto(rv_unit, out);
        } else if (comptime U_in.dim.eql(frequency_dim)) {
            var temp_f: Quantity(T, si.Hz) = .init(out.value);
            input_quantity.toInto(si.Hz, &temp_f);
            fToRvInto(T, &temp_f, rest, out);
        } else if (comptime U_in.dim.eql(wavelength_dim) or
            U_in.dim.eql(energy_dim) or
            U_in.dim.eql(wavenumber_spec_dim) or
            U_in.dim.eql(wavenumber_ang_dim))
        {
            var temp_f: Quantity(T, si.Hz) = .init(out.value);
            spectral_equivalency.convertInto(T, U_in, input_quantity, si.Hz, &temp_f);
            fToRvInto(T, &temp_f, rest, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromRv(comptime T: type, input_quantity: *const Quantity(T, rv_unit), comptime U_out: Unit, rest: *const Quantity(f64, si.Hz)) Quantity(T, U_out) {
        if (comptime U_out.dim.eql(rv_dim)) {
            return input_quantity.to(U_out);
        } else if (comptime U_out.dim.eql(frequency_dim)) {
            const f = rvToF(T, input_quantity, rest);
            return f.to(U_out);
        } else if (comptime U_out.dim.eql(energy_dim) or
            U_out.dim.eql(wavelength_dim) or
            U_out.dim.eql(wavenumber_spec_dim) or
            U_out.dim.eql(wavenumber_ang_dim))
        {
            const f = rvToF(T, input_quantity, rest);
            return spectral_equivalency.convert(T, si.Hz, &f, U_out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromRvInto(comptime T: type, input_quantity: *const Quantity(T, rv_unit), comptime U_out: Unit, rest: *const Quantity(f64, si.Hz), out: *Quantity(T, U_out)) void {
        if (comptime U_out.dim.eql(rv_dim)) {
            input_quantity.toInto(U_out, out);
        } else if (comptime U_out.dim.eql(frequency_dim)) {
            var temp_f: Quantity(T, si.Hz) = .init(out.value);
            rvToFInto(T, input_quantity, rest, &temp_f);
            temp_f.toInto(U_out, out);
        } else if (comptime U_out.dim.eql(energy_dim) or
            U_out.dim.eql(wavelength_dim) or
            U_out.dim.eql(wavenumber_spec_dim) or
            U_out.dim.eql(wavenumber_ang_dim))
        {
            var temp_f: Quantity(T, si.Hz) = .init(out.value);
            rvToFInto(T, input_quantity, rest, &temp_f);
            spectral_equivalency.convertInto(T, si.Hz, &temp_f, U_out, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fToRv(comptime T: type, freq: *const Quantity(T, si.Hz), rest: *const Quantity(f64, si.Hz)) Quantity(T, rv_unit) {
        // v = c * (f0 - f) / f0
        return (freq.neg().addScalar(rest.*)).divScalar(rest.*).mulScalar(c_kms).to(rv_unit);
    }

    fn fToRvInto(comptime T: type, freq: *const Quantity(T, si.Hz), rest: *const Quantity(f64, si.Hz), out: *Quantity(T, rv_unit)) void {
        const f0 = rest.*;
        const inv_f0 = f0.pow(-1);

        // V/c = 1 - f/f_0
        var ratio: Quantity(T, @TypeOf(freq.*).unit.mul(@TypeOf(inv_f0).unit)) = .init(out.value);
        freq.mulScalarInto(inv_f0, &ratio); // f/f0
        var ratio_unitless: Quantity([]f64, UNITLESS) = .init(out.value);
        ratio.toInto(UNITLESS, &ratio_unitless);
        ratio_unitless.mulScalarValueInPlace(-1); // -f/f0
        ratio_unitless.addScalarValueInPlace(1); // 1 - f/f0 = (f_0 - f)/f_0
        var vel: Quantity(T, UNITLESS.mul(rv_unit)) = .init(out.value);
        ratio_unitless.mulScalarInto(c_kms, &vel); // c . (f_0 - f)/f_0 = V
        vel.toInto(rv_unit, out);
    }

    fn rvToF(comptime T: type, rv: *const Quantity(T, rv_unit), rest: *const Quantity(f64, si.Hz)) Quantity(T, si.Hz) {
        // f = f0 * (1 - v/c)
        const beta = rv.divScalar(c_kms).to(UNITLESS);
        return (beta.neg().addScalarValue(1)).mulScalar(rest.*);
    }

    fn rvToFInto(comptime T: type, rv: *const Quantity(T, rv_unit), rest: *const Quantity(f64, si.Hz), out: *Quantity(T, si.Hz)) void {
        const f0 = rest.*;
        // f = f0 . (1 - beta), beta = V/c
        var beta_raw: Quantity(T, rv_unit.div(rv_unit)) = .init(out.value);
        rv.divScalarInto(c_kms, &beta_raw); // V/c
        var buf: Quantity(T, UNITLESS) = .init(out.value);
        beta_raw.toInto(UNITLESS, &buf); // beta
        buf.mulScalarValueInPlace(-1); // -beta
        buf.addScalarValueInPlace(1); // 1 - beta
        var freq: Quantity(T, UNITLESS.mul(si.Hz)) = .init(out.value);
        buf.mulScalarInto(f0, &freq); // f_0 . (1 - beta)
        freq.toInto(si.Hz, out);
    }
};

/// Magnetic flux-field equivalency: converts between the H-field (A/m) and the
/// B-field (tesla) via B = mu_0 mu_r H. Takes an optional argument: the relative
/// permeability mu_r as a dimensionless Quantity (defaults to 1, i.e. vacuum).
/// Hub unit is the tesla.
pub const magnetic_flux_field_equivalency = struct {
    const name = "magnetic-flux";
    const H_field_dim: Dim = si.A.div(si.m).dim;
    const B_field_dim: Dim = si.T.dim;
    const mu_0 = cst.mu0.quantity;

    fn checkArgsAndGetMuR(args: anytype) Quantity(f64, UNITLESS) {
        const ArgsType = @TypeOf(args);
        const args_fields = @typeInfo(ArgsType).@"struct".fields;

        if (comptime args_fields.len > 1) {
            @compileError("Magnetic-flux field equivalency expect at most 1 argument.");
        }

        if (comptime args_fields.len == 0) {
            return .init(1);
        }

        const ArgType = args_fields[0].type;
        utils.assertIsScalarQuantity(ArgType);
        if (comptime !ArgType.unit.dim.eql(UNITLESS.dim)) {
            @compileError("Magnetic-flux field equivalency argument must be a unitless Quantity.");
        }
        const value_f64: f64 = @floatCast(args.@"0".value);
        return .init(value_f64);
    }

    pub fn convert(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, args: anytype) Quantity(T, U_out) {
        const mu_r = checkArgsAndGetMuR(args);
        const as_T = toT(T, U_in, input_quantity, &mu_r);
        return fromT(T, &as_T, U_out, &mu_r);
    }

    pub fn convertInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, out: *Quantity(T, U_out), args: anytype) void {
        const mu_r = checkArgsAndGetMuR(args);
        var temp_T: Quantity(T, si.T) = .init(out.value);
        toTInto(T, U_in, input_quantity, &mu_r, &temp_T);
        fromTInto(T, &temp_T, U_out, &mu_r, out);
    }

    fn toT(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), mu_r: *const Quantity(f64, UNITLESS)) Quantity(T, si.T) {
        if (comptime U_in.dim.eql(B_field_dim)) {
            return input_quantity.to(si.T);
        } else if (comptime U_in.dim.eql(H_field_dim)) {
            return input_quantity.mulScalar(mu_0).mulScalar(mu_r.*).to(si.T);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn toTInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), mu_r: *const Quantity(f64, UNITLESS), out: *Quantity(T, si.T)) void {
        if (comptime U_in.dim.eql(B_field_dim)) {
            input_quantity.toInto(si.T, out);
        } else if (comptime U_in.dim.eql(H_field_dim)) {
            const mu_0_mu_r = mu_0.mulScalar(mu_r.*);
            var temp: Quantity(T, U_in.mul(@TypeOf(mu_0_mu_r).unit)) = .init(out.value);
            input_quantity.mulScalarInto(mu_0_mu_r, &temp);
            temp.toInto(si.T, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromT(comptime T: type, input_quantity: *const Quantity(T, si.T), comptime U_out: Unit, mu_r: *const Quantity(f64, UNITLESS)) Quantity(T, U_out) {
        if (comptime U_out.dim.eql(B_field_dim)) {
            return input_quantity.to(U_out);
        } else if (comptime U_out.dim.eql(H_field_dim)) {
            return input_quantity.divScalar(mu_0.mulScalar(mu_r.*)).to(U_out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromTInto(comptime T: type, input_quantity: *const Quantity(T, si.T), comptime U_out: Unit, mu_r: *const Quantity(f64, UNITLESS), out: *Quantity(T, U_out)) void {
        if (comptime U_out.dim.eql(B_field_dim)) {
            input_quantity.toInto(U_out, out);
        } else if (comptime U_out.dim.eql(H_field_dim)) {
            const mu_0_mu_r = mu_0.mulScalar(mu_r.*);
            var temp: Quantity(T, si.T.div(@TypeOf(mu_0_mu_r).unit)) = .init(out.value);
            input_quantity.divScalarInto(mu_0_mu_r, &temp);
            temp.toInto(U_out, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }
};

// NOTE: In astropy, this conversion also support Jy/beam unit. In Zali at the moment we do not
// support beam unit which is a dimensionless unit that seems misleading to me.
/// Brightness-temperature equivalency: converts between specific intensity / flux
/// density and brightness temperature using the Rayleigh-Jeans law. Hub unit is the
/// kelvin. Takes one or two arguments: the spectral location (frequency, required),
/// and a beam solid angle (required only when the flux side is per-beam Jy rather
/// than per-steradian Jy/sr).
///   S_nu = 2 k_B T_B nu^2 Omega / c^2   (Omega = beam solid angle)
pub const brightness_temperature_equivalency = struct {
    const name = "brightness-temperature";
    const sterad_dim: Dim = si.sr.dim;
    const Jy: Unit = us.Jy;
    const Jy_dim: Dim = Jy.dim;
    const K_dim: Dim = si.K.dim;
    const JySr: Unit = us.Jy.div(si.sr);
    const JySr_dim: Dim = JySr.dim;
    const c2 = cst.c.quantity.pow(2);
    const two_kb = cst.k_B.quantity.mulScalarValue(2);
    const one_sr: Quantity(f64, si.sr) = .init(1);

    fn checkArgs(args: anytype) void {
        const ArgsType = @TypeOf(args);
        const args_fields = @typeInfo(ArgsType).@"struct".fields;

        if (comptime args_fields.len != 1 and args_fields.len != 2) {
            @compileError("Brightness-temperature equivalency expect 1 or 2 argument.");
        }
        inline for (args_fields) |field| {
            const ArgType = field.type;
            utils.assertIsScalarQuantity(ArgType);
        }

        if (comptime !utils.isSpectralUnit(args_fields[0].type.unit)) {
            @compileError("Brightness-temperature equivalency first argument must be a spectral compatible Quantity.");
        }
        if (comptime args_fields.len == 2 and !args_fields[1].type.unit.dim.eql(sterad_dim)) {
            @compileError("Brightness-temperature equivalency second arguement must be in square angular units, i.e. steradian equivalent.");
        }
    }

    pub fn convert(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, args: anytype) Quantity(T, U_out) {
        checkArgs(args);
        const arg_number: usize = @typeInfo(@TypeOf(args)).@"struct".fields.len;
        const spectral_quantity = args.@"0";
        const spectral_quantity_unit = comptime @TypeOf(spectral_quantity).unit;
        const spectral_quantity_value_type = comptime @TypeOf(spectral_quantity).value_type;
        const freq: Quantity(spectral_quantity_value_type, si.Hz) = spectral_equivalency.convert(spectral_quantity_value_type, spectral_quantity_unit, &spectral_quantity, si.Hz);
        const freq_value_f64: f64 = @floatCast(freq.value);
        const freq_f64: Quantity(f64, si.Hz) = .init(freq_value_f64);
        if (comptime arg_number == 2) {
            const beam = args.@"1".to(si.sr);
            const beam_value_f64: f64 = @floatCast(beam.value);
            const beam_f64: Quantity(f64, si.sr) = .init(beam_value_f64);
            const as_K = toK(T, U_in, arg_number, input_quantity, &freq_f64, &beam_f64);
            return fromK(T, &as_K, U_out, arg_number, &freq_f64, &beam_f64);
        } else {
            const as_K = toK(T, U_in, arg_number, input_quantity, &freq_f64, {});
            return fromK(T, &as_K, U_out, arg_number, &freq_f64, {});
        }
    }

    pub fn convertInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, out: *Quantity(T, U_out), args: anytype) void {
        checkArgs(args);
        const arg_number: usize = @typeInfo(@TypeOf(args)).@"struct".fields.len;
        const spectral_quantity = args.@"0";
        const spectral_quantity_unit = comptime @TypeOf(spectral_quantity).unit;
        const spectral_quantity_value_type = comptime @TypeOf(spectral_quantity).value_type;
        const freq: Quantity(spectral_quantity_value_type, si.Hz) = spectral_equivalency.convert(spectral_quantity_value_type, spectral_quantity_unit, &spectral_quantity, si.Hz);
        const freq_value_f64: f64 = @floatCast(freq.value);
        const freq_f64: Quantity(f64, si.Hz) = .init(freq_value_f64);
        var temp_K: Quantity(T, si.K) = .init(out.value);
        if (comptime arg_number == 2) {
            const beam = args.@"1".to(si.sr);
            const beam_value_f64: f64 = @floatCast(beam.value);
            const beam_f64: Quantity(f64, si.sr) = .init(beam_value_f64);
            toKInto(T, U_in, arg_number, input_quantity, &freq_f64, &beam_f64, &temp_K);
            fromKInto(T, &temp_K, U_out, arg_number, &freq_f64, &beam_f64, out);
        } else {
            toKInto(T, U_in, arg_number, input_quantity, &freq_f64, {}, &temp_K);
            fromKInto(T, &temp_K, U_out, arg_number, &freq_f64, {}, out);
        }
    }

    fn toK(
        comptime T: type,
        comptime U_in: Unit,
        comptime n_args: usize,
        input_quantity: *const Quantity(T, U_in),
        freq: *const Quantity(f64, si.Hz),
        beam: if (n_args == 2) *const Quantity(f64, si.sr) else void,
    ) Quantity(T, si.K) {
        const nu2 = freq.pow(2);
        if (comptime U_in.dim.eql(K_dim)) {
            return input_quantity.to(si.K);
        } else if (comptime U_in.dim.eql(Jy_dim)) {
            if (comptime n_args != 2) @compileError("Converting Jy to K required a beam area as a second argument.");
            return input_quantity.to(Jy).mulScalar(c2).divScalar(two_kb.mul(nu2).mul(beam.*)).mulScalar(one_sr).to(si.K);
        } else if (comptime U_in.dim.eql(JySr_dim)) {
            return input_quantity.to(JySr).mulScalar(c2).div(two_kb.mul(nu2)).mulScalar(one_sr).to(si.K);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn toKInto(
        comptime T: type,
        comptime U_in: Unit,
        comptime n_args: usize,
        input_quantity: *const Quantity(T, U_in),
        freq: *const Quantity(f64, si.Hz),
        beam: if (n_args == 2) *const Quantity(f64, si.sr) else void,
        out: *Quantity(T, si.K),
    ) void {
        const nu2 = freq.pow(2);
        if (comptime U_in.dim.eql(K_dim)) {
            input_quantity.toInto(si.K, out);
        } else if (comptime U_in.dim.eql(Jy_dim)) {
            if (comptime n_args != 2) @compileError("Converting Jy to K required a beam area as a second argument.");
            var temp_Jy: Quantity(T, Jy) = .init(out.value);
            input_quantity.toInto(Jy, &temp_Jy);
            const scalar = c2.div(two_kb.mul(nu2).mul(beam.*)).mul(one_sr);
            var temp_temp: Quantity(T, Jy.mul(@TypeOf(scalar).unit)) = .init(out.value);
            temp_Jy.mulScalarInto(scalar, &temp_temp);
            temp_temp.toInto(si.K, out);
        } else if (comptime U_in.dim.eql(JySr_dim)) {
            var temp_JySr: Quantity(T, JySr) = .init(out.value);
            input_quantity.toInto(JySr, &temp_JySr);
            const scalar = c2.div(two_kb.mul(nu2)).mul(one_sr);
            var temp_temp: Quantity(T, JySr.mul(@TypeOf(scalar).unit)) = .init(out.value);
            temp_JySr.mulScalarInto(scalar, &temp_temp);
            temp_temp.toInto(si.K, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromK(
        comptime T: type,
        input_quantity: *const Quantity(T, si.K),
        comptime U_out: Unit,
        comptime n_args: usize,
        freq: *const Quantity(f64, si.Hz),
        beam: if (n_args == 2) *const Quantity(f64, si.sr) else void,
    ) Quantity(T, U_out) {
        const nu2 = freq.pow(2);
        if (comptime U_out.dim.eql(K_dim)) {
            return input_quantity.to(U_out);
        } else if (comptime U_out.dim.eql(Jy_dim)) {
            if (comptime n_args != 2) @compileError("Converting K to Jy required a beam area as a second argument.");
            return input_quantity.mulScalar(two_kb.mul(nu2).mul(beam.*).div(c2).div(one_sr)).to(U_out);
        } else if (comptime U_out.dim.eql(JySr_dim)) {
            return input_quantity.mulScalar(two_kb.mul(nu2).div(c2).div(one_sr)).to(U_out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromKInto(
        comptime T: type,
        input_quantity: *const Quantity(T, si.K),
        comptime U_out: Unit,
        comptime n_args: usize,
        freq: *const Quantity(f64, si.Hz),
        beam: if (n_args == 2) *const Quantity(f64, si.sr) else void,
        out: *Quantity(T, U_out),
    ) void {
        const nu2 = freq.pow(2);
        if (comptime U_out.dim.eql(K_dim)) {
            input_quantity.toInto(U_out, out);
        } else if (comptime U_out.dim.eql(Jy_dim)) {
            if (comptime n_args != 2) @compileError("Converting K to Jy required a beam area as a second argument.");
            const scalar = two_kb.mul(nu2).mul(beam.*).div(c2).div(one_sr);
            var temp: Quantity(T, si.K.mul(@TypeOf(scalar).unit)) = .init(out.value);
            input_quantity.mulScalarInto(scalar, &temp);
            temp.toInto(U_out, out);
        } else if (comptime U_out.dim.eql(JySr_dim)) {
            const scalar = two_kb.mul(nu2).div(c2).div(one_sr);
            var temp: Quantity(T, si.K.mul(@TypeOf(scalar).unit)) = .init(out.value);
            input_quantity.mulScalarInto(scalar, &temp);
            temp.toInto(U_out, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }
};

/// Plate-scale equivalency: converts between a focal-plane length and an angle on the
/// sky, given a plate scale in length/angle (e.g. m/rad). Direct two-way conversion.
///   length = angle * plate_scale,   angle = length / plate_scale
pub const plate_scale_equivalency = struct {
    const name = "plate-scale";
    const plate_scale_dim: Dim = si.m.div(si.rad).dim;
    const length_dim: Dim = dim_mod.length;
    const angle_dim: Dim = dim_mod.angle;

    fn checkArgs(args: anytype) void {
        const ArgsType = @TypeOf(args);
        const args_fields = @typeInfo(ArgsType).@"struct".fields;

        if (comptime args_fields.len != 1) {
            @compileError("Plate-scale equivalency requires one argument.");
        }

        const ArgType = args_fields[0].type;
        utils.assertIsScalarQuantity(ArgType);
        if (comptime !ArgType.unit.dim.eql(plate_scale_dim)) {
            @compileError("Plate-scale equivalency argument must be in units of length / angle.");
        }
    }

    pub fn convert(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, args: anytype) Quantity(T, U_out) {
        checkArgs(args);
        const plate_scale = args.@"0";
        const plate_scale_unit = comptime @TypeOf(plate_scale).unit;
        const plate_scale_value_f64: f64 = @floatCast(plate_scale.value);
        const plate_scale_f64: Quantity(f64, plate_scale_unit) = .init(plate_scale_value_f64);
        if (comptime U_in.dim.eql(length_dim) and U_out.dim.eql(angle_dim)) {
            return input_quantity.divScalar(plate_scale_f64).to(U_out);
        } else if (comptime U_in.dim.eql(angle_dim) and U_out.dim.eql(length_dim)) {
            return input_quantity.mulScalar(plate_scale_f64).to(U_out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    pub fn convertInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, out: *Quantity(T, U_out), args: anytype) void {
        checkArgs(args);
        const plate_scale = args.@"0";
        const plate_scale_unit = comptime @TypeOf(plate_scale).unit;
        const plate_scale_value_f64: f64 = @floatCast(plate_scale.value);
        const plate_scale_f64: Quantity(f64, plate_scale_unit) = .init(plate_scale_value_f64);
        if (comptime U_in.dim.eql(length_dim) and U_out.dim.eql(angle_dim)) {
            var temp: Quantity(T, U_in.div(plate_scale_unit)) = .init(out.value);
            input_quantity.divScalarInto(plate_scale_f64, &temp);
            temp.toInto(U_out, out);
        } else if (comptime U_in.dim.eql(angle_dim) and U_out.dim.eql(length_dim)) {
            var temp: Quantity(T, U_in.mul(plate_scale_unit)) = .init(out.value);
            input_quantity.mulScalarInto(plate_scale_f64, &temp);
            temp.toInto(U_out, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }
};

/// Pixel-scale equivalency: converts between pixels and a physical quantity (an angle
/// or a length), given a pixel scale such as arcsec/pixel or pixel/cm. Because pixel
/// is dimensionless, the conversion is driven by inspecting the pixel exponent in the
/// scale's symbol rather than by dimension alone.
pub const pixel_scale_equivalency = struct {
    const name = "pixel-scale";

    fn checkArgs(args: anytype) void {
        const ArgsType = @TypeOf(args);
        const args_fields = @typeInfo(ArgsType).@"struct".fields;

        if (comptime args_fields.len != 1) {
            @compileError("Pixel-scale equivalency requires one argument.");
        }

        const ArgType = args_fields[0].type;
        utils.assertIsScalarQuantity(ArgType);
    }

    fn getPixelPower(comptime unit: Unit) isize {
        const symbol = unit.symbol;
        if (comptime symbol.indexOfSymbol(misc.pix.symbol.terms[0])) |idx| {
            const pixel_power_frac = symbol.terms[idx].exponent;
            if (comptime pixel_power_frac.eqlScalar(-1)) return -1;
            if (comptime pixel_power_frac.eqlScalar(1)) return 1;
            @compileError("Pixel-scale equivalency argument must have pixel dimensionality of 1 or -1.");
        } else {
            @compileError("Pixel-scale equivalency argument must contain a pixel unit.");
        }
    }

    pub fn convert(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, args: anytype) Quantity(T, U_out) {
        checkArgs(args);
        const pixel_scale = args.@"0";
        const pixel_scale_type = comptime @TypeOf(pixel_scale);
        const pixel_scale_unit = pixel_scale_type.unit;
        const pixel_scale_value_f64: f64 = @floatCast(pixel_scale.value);
        const pixel_scale_f64: Quantity(f64, pixel_scale_unit) = .init(pixel_scale_value_f64);
        const pixel_power = comptime getPixelPower(pixel_scale_unit);

        const physical_unit: Unit = comptime blk: {
            if (pixel_power == -1) break :blk pixel_scale_unit.mul(misc.pix);
            break :blk pixel_scale_unit.div(misc.pix).pow(-1);
        };

        const scale: Quantity(f64, UNITLESS) = .initScalarValue(if (comptime pixel_power == -1)
            pixel_scale_f64.to(physical_unit.div(misc.pix)).value
        else
            pixel_scale_f64.to(misc.pix.div(physical_unit)).value);

        // pix → physical
        if (comptime U_in.dim.eql(misc.pix.dim) and U_out.dim.eql(physical_unit.dim)) {
            if (comptime pixel_power == -1) {
                return Quantity(T, physical_unit).init(input_quantity.value).mulScalar(scale).to(U_out);
            } else {
                return Quantity(T, physical_unit).init(input_quantity.value).divScalar(scale).to(U_out);
            }
        }
        // physical → pix
        else if (comptime U_in.dim.eql(physical_unit.dim) and U_out.dim.eql(misc.pix.dim)) {
            if (comptime pixel_power == -1) {
                return Quantity(T, misc.pix).init(input_quantity.to(physical_unit).value).divScalar(scale).to(U_out);
            } else {
                return Quantity(T, misc.pix).init(input_quantity.to(physical_unit).value).mulScalar(scale).to(U_out);
            }
        } else {
            @compileError(unknownConversion(name));
        }
    }

    pub fn convertInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, out: *Quantity(T, U_out), args: anytype) void {
        checkArgs(args);
        const pixel_scale = args.@"0";
        const pixel_scale_type = comptime @TypeOf(pixel_scale);
        const pixel_scale_unit = pixel_scale_type.unit;
        const pixel_scale_value_f64: f64 = @floatCast(pixel_scale.value);
        const pixel_scale_f64: Quantity(f64, pixel_scale_unit) = .init(pixel_scale_value_f64);
        const pixel_power = comptime getPixelPower(pixel_scale_unit);

        const physical_unit: Unit = comptime blk: {
            if (pixel_power == -1) break :blk pixel_scale_unit.mul(misc.pix);
            break :blk pixel_scale_unit.div(misc.pix).pow(-1);
        };

        const scale: Quantity(f64, UNITLESS) = .initScalarValue(if (comptime pixel_power == -1)
            pixel_scale_f64.to(physical_unit.div(misc.pix)).value
        else
            pixel_scale_f64.to(misc.pix.div(physical_unit)).value);

        // pix → physical
        if (comptime U_in.dim.eql(misc.pix.dim) and U_out.dim.eql(physical_unit.dim)) {
            if (comptime pixel_power == -1) {
                var temp: Quantity(T, physical_unit.mul(UNITLESS)) = .init(out.value);
                @memcpy(out.value, input_quantity.value);
                temp.mulScalarInto(scale, &temp);
                temp.toInto(U_out, out);
            } else {
                var temp: Quantity(T, physical_unit.div(UNITLESS)) = .init(out.value);
                @memcpy(out.value, input_quantity.value);
                temp.divScalarInto(scale, &temp);
                temp.toInto(U_out, out);
            }
        }
        // physical → pix
        else if (comptime U_in.dim.eql(physical_unit.dim) and U_out.dim.eql(misc.pix.dim)) {
            if (comptime pixel_power == -1) {
                var physical: Quantity(T, physical_unit) = .init(out.value);
                input_quantity.toInto(physical_unit, &physical);
                var temp: Quantity(T, misc.pix) = .init(out.value);
                temp.divScalarInto(scale, &temp);
                temp.toInto(U_out, out);
            } else {
                var physical: Quantity(T, physical_unit) = .init(out.value);
                input_quantity.toInto(physical_unit, &physical);
                var temp: Quantity(T, misc.pix) = .init(out.value);
                temp.mulScalarInto(scale, &temp);
                temp.toInto(U_out, out);
            }
        } else {
            @compileError(unknownConversion(name));
        }
    }
};

/// Spectral-density equivalency: converts between per-frequency and per-wavelength
/// spectral densities (and their nu*F_nu integrated form) at a given spectral
/// location, applying the Jacobian dnu/dlambda = c/lambda^2. Supports four families,
/// each with its own per-frequency hub:
///   - flux density           (F_nu  <-> F_lambda  <-> nu F_nu)
///   - luminosity density     (L_nu  <-> L_lambda  <-> nu L_nu)
///   - surface brightness     (S_nu  <-> S_lambda  <-> nu S_nu),      per steradian
///   - luminosity surface br. (SL_nu <-> SL_lambda <-> nu SL_nu),     per steradian
/// The Jacobian is identical across families; only the extra dimensions differ.
/// Requires one argument: the spectral location as a spectral-compatible Quantity.
// TODO: support photon base units
pub const spectral_density_equivalency = struct {
    const name = "spectral-density";
    const c = cst.c.quantity;
    // Flux
    const Fnu: Unit = si.kg.div(si.s.pow(2));
    const Flamb: Unit = si.kg.div(si.m.mul(si.s.pow(3)));
    const nuFnu: Unit = si.kg.div(si.s.pow(3));
    const Fnu_dim: Dim = Fnu.dim;
    const Flamb_dim: Dim = Flamb.dim;
    const nuFnu_dim: Dim = nuFnu.dim;

    // Luminosity
    const Lnu: Unit = si.W.div(si.Hz);
    const Llamb: Unit = si.W.div(si.m);
    const nuLnu: Unit = si.W;
    const Lnu_dim: Dim = Lnu.dim;
    const Llamb_dim: Dim = Llamb.dim;
    const nuLnu_dim: Dim = nuLnu.dim;

    // Surface brightness
    const Snu: Unit = si.W.div(si.m.pow(2)).div(si.Hz).div(si.sr);
    const Slamb: Unit = si.W.div(si.m.pow(3)).div(si.sr);
    const nuSnu: Unit = si.W.div(si.m.pow(2)).div(si.sr);
    const Snu_dim: Dim = Snu.dim;
    const Slamb_dim: Dim = Slamb.dim;
    const nuSnu_dim: Dim = nuSnu.dim;

    // Luminosity surface brightness
    const SLnu: Unit = si.W.div(si.Hz).div(si.sr);
    const SLlamb: Unit = si.W.div(si.m).div(si.sr);
    const nuSLnu: Unit = si.W.div(si.sr);
    const SLnu_dim: Dim = SLnu.dim;
    const SLlamb_dim: Dim = SLlamb.dim;
    const nuSLnu_dim: Dim = nuSLnu.dim;

    fn checkArgs(args: anytype) void {
        const ArgsType = @TypeOf(args);
        const args_fields = @typeInfo(ArgsType).@"struct".fields;

        if (comptime args_fields.len != 1) {
            @compileError("Spectral-density equivalency expect one argument.");
        }

        const ArgType = args_fields[0].type;
        utils.assertIsScalarQuantity(ArgType);
        if (comptime !utils.isSpectralUnit(ArgType.unit)) {
            @compileError("Spectral-density equivalency argument must be a spectral compatible Quantity.");
        }
    }

    pub fn convert(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, args: anytype) Quantity(T, U_out) {
        checkArgs(args);
        const spectral_quantity = args.@"0";
        const spectral_quantity_unit = comptime @TypeOf(spectral_quantity).unit;
        const spectral_quantity_value_f64: f64 = @floatCast(spectral_quantity.value);
        const spectral_quantity_f64: Quantity(f64, spectral_quantity_unit) = .init(spectral_quantity_value_f64);
        const freq: Quantity(f64, si.Hz) = spectral_equivalency.convert(f64, spectral_quantity_unit, &spectral_quantity_f64, si.Hz);
        if (comptime U_in.dim.eql(Fnu_dim) or U_in.dim.eql(Flamb_dim) or U_in.dim.eql(nuFnu_dim)) {
            const as_Fnu = toFnu(T, U_in, input_quantity, &freq);
            return fromFnu(T, &as_Fnu, U_out, &freq);
        } else if (comptime U_in.dim.eql(Lnu_dim) or U_in.dim.eql(Llamb_dim) or U_in.dim.eql(nuLnu_dim)) {
            const as_Lnu = toLnu(T, U_in, input_quantity, &freq);
            return fromLnu(T, &as_Lnu, U_out, &freq);
        } else if (comptime U_in.dim.eql(Snu_dim) or U_in.dim.eql(Slamb_dim) or U_in.dim.eql(nuSnu_dim)) {
            const as_Snu = toSnu(T, U_in, input_quantity, &freq);
            return fromSnu(T, &as_Snu, U_out, &freq);
        } else if (comptime U_in.dim.eql(SLnu_dim) or U_in.dim.eql(SLlamb_dim) or U_in.dim.eql(nuSLnu_dim)) {
            const as_SLnu = toSLnu(T, U_in, input_quantity, &freq);
            return fromSLnu(T, &as_SLnu, U_out, &freq);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    pub fn convertInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, out: *Quantity(T, U_out), args: anytype) void {
        checkArgs(args);
        const spectral_quantity = args.@"0";
        const spectral_quantity_unit = comptime @TypeOf(spectral_quantity).unit;
        const spectral_quantity_value_f64: f64 = @floatCast(spectral_quantity.value);
        const spectral_quantity_f64: Quantity(f64, spectral_quantity_unit) = .init(spectral_quantity_value_f64);
        const freq: Quantity(f64, si.Hz) = spectral_equivalency.convert(f64, spectral_quantity_unit, &spectral_quantity_f64, si.Hz);
        if (comptime U_in.dim.eql(Fnu_dim) or U_in.dim.eql(Flamb_dim) or U_in.dim.eql(nuFnu_dim)) {
            var temp: Quantity(T, Fnu) = .init(out.value);
            toFnuInto(T, U_in, input_quantity, &freq, &temp);
            fromFnuInto(T, &temp, U_out, &freq, out);
        } else if (comptime U_in.dim.eql(Lnu_dim) or U_in.dim.eql(Llamb_dim) or U_in.dim.eql(nuLnu_dim)) {
            var temp: Quantity(T, Lnu) = .init(out.value);
            toLnuInto(T, U_in, input_quantity, &freq, &temp);
            fromLnuInto(T, &temp, U_out, &freq, out);
        } else if (comptime U_in.dim.eql(Snu_dim) or U_in.dim.eql(Slamb_dim) or U_in.dim.eql(nuSnu_dim)) {
            var temp: Quantity(T, Snu) = .init(out.value);
            toSnuInto(T, U_in, input_quantity, &freq, &temp);
            fromSnuInto(T, &temp, U_out, &freq, out);
        } else if (comptime U_in.dim.eql(SLnu_dim) or U_in.dim.eql(SLlamb_dim) or U_in.dim.eql(nuSLnu_dim)) {
            var temp: Quantity(T, SLnu) = .init(out.value);
            toSLnuInto(T, U_in, input_quantity, &freq, &temp);
            fromSLnuInto(T, &temp, U_out, &freq, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn toFnu(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), freq: *const Quantity(f64, si.Hz)) Quantity(T, Fnu) {
        if (comptime U_in.dim.eql(Fnu_dim)) {
            return input_quantity.to(Fnu);
        } else if (comptime U_in.dim.eql(Flamb_dim)) {
            return input_quantity.mulScalar(c).divScalar(freq.pow(2)).to(Fnu);
        } else if (comptime U_in.dim.eql(nuFnu_dim)) {
            return input_quantity.divScalar(freq.*).to(Fnu);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn toFnuInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), freq: *const Quantity(f64, si.Hz), out: *Quantity(T, Fnu)) void {
        if (comptime U_in.dim.eql(Fnu_dim)) {
            input_quantity.toInto(Fnu, out);
        } else if (comptime U_in.dim.eql(Flamb_dim)) {
            const scalar = c.div(freq.pow(2));
            var temp: Quantity(T, U_in.mul(@TypeOf(scalar).unit)) = .init(out.value);
            input_quantity.mulScalarInto(scalar, &temp);
            temp.toInto(Fnu, out);
        } else if (comptime U_in.dim.eql(nuFnu_dim)) {
            var temp: Quantity(T, U_in.div(@TypeOf(freq.*).unit)) = .init(out.value);
            input_quantity.divScalarInto(freq.*, &temp);
            temp.toInto(Fnu, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromFnu(comptime T: type, input_quantity: *const Quantity(T, Fnu), comptime U_out: Unit, freq: *const Quantity(f64, si.Hz)) Quantity(T, U_out) {
        if (comptime U_out.dim.eql(Fnu_dim)) {
            return input_quantity.to(U_out);
        } else if (comptime U_out.dim.eql(Flamb_dim)) {
            return input_quantity.mulScalar(freq.pow(2)).divScalar(c).to(U_out);
        } else if (comptime U_out.dim.eql(nuFnu_dim)) {
            return input_quantity.mulScalar(freq.*).to(U_out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromFnuInto(comptime T: type, input_quantity: *const Quantity(T, Fnu), comptime U_out: Unit, freq: *const Quantity(f64, si.Hz), out: *Quantity(T, U_out)) void {
        if (comptime U_out.dim.eql(Fnu_dim)) {
            input_quantity.toInto(U_out, out);
        } else if (comptime U_out.dim.eql(Flamb_dim)) {
            const scalar = freq.pow(2).div(c);
            var temp: Quantity(T, Fnu.mul(@TypeOf(scalar).unit)) = .init(out.value);
            input_quantity.mulScalarInto(scalar, &temp);
            temp.toInto(U_out, out);
        } else if (comptime U_out.dim.eql(nuFnu_dim)) {
            var temp: Quantity(T, Fnu.mul(@TypeOf(freq.*).unit)) = .init(out.value);
            input_quantity.mulScalarInto(freq.*, &temp);
            temp.toInto(U_out, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn toLnu(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), freq: *const Quantity(f64, si.Hz)) Quantity(T, Lnu) {
        if (comptime U_in.dim.eql(Lnu_dim)) {
            return input_quantity.to(Lnu);
        } else if (comptime U_in.dim.eql(Llamb_dim)) {
            return input_quantity.mulScalar(c).divScalar(freq.pow(2)).to(Lnu);
        } else if (comptime U_in.dim.eql(nuLnu_dim)) {
            return input_quantity.divScalar(freq.*).to(Lnu);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn toLnuInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), freq: *const Quantity(f64, si.Hz), out: *Quantity(T, Lnu)) void {
        if (comptime U_in.dim.eql(Lnu_dim)) {
            input_quantity.toInto(Lnu, out);
        } else if (comptime U_in.dim.eql(Llamb_dim)) {
            const scalar = c.div(freq.pow(2));
            var temp: Quantity(T, U_in.mul(@TypeOf(scalar).unit)) = .init(out.value);
            input_quantity.mulScalarInto(scalar, &temp);
            temp.toInto(Lnu, out);
        } else if (comptime U_in.dim.eql(nuLnu_dim)) {
            var temp: Quantity(T, U_in.div(@TypeOf(freq.*).unit)) = .init(out.value);
            input_quantity.divScalarInto(freq.*, &temp);
            temp.toInto(Lnu, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromLnu(comptime T: type, input_quantity: *const Quantity(T, Lnu), comptime U_out: Unit, freq: *const Quantity(f64, si.Hz)) Quantity(T, U_out) {
        if (comptime U_out.dim.eql(Lnu_dim)) {
            return input_quantity.to(U_out);
        } else if (comptime U_out.dim.eql(Llamb_dim)) {
            return input_quantity.mulScalar(freq.pow(2)).divScalar(c).to(U_out);
        } else if (comptime U_out.dim.eql(nuLnu_dim)) {
            return input_quantity.mulScalar(freq.*).to(U_out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromLnuInto(comptime T: type, input_quantity: *const Quantity(T, Lnu), comptime U_out: Unit, freq: *const Quantity(f64, si.Hz), out: *Quantity(T, U_out)) void {
        if (comptime U_out.dim.eql(Lnu_dim)) {
            input_quantity.toInto(U_out, out);
        } else if (comptime U_out.dim.eql(Llamb_dim)) {
            const scalar = freq.pow(2).div(c);
            var temp: Quantity(T, Lnu.mul(@TypeOf(scalar).unit)) = .init(out.value);
            input_quantity.mulScalarInto(scalar, &temp);
            temp.toInto(U_out, out);
        } else if (comptime U_out.dim.eql(nuLnu_dim)) {
            var temp: Quantity(T, Lnu.mul(@TypeOf(freq.*).unit)) = .init(out.value);
            input_quantity.mulScalarInto(freq.*, &temp);
            temp.toInto(U_out, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn toSnu(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), freq: *const Quantity(f64, si.Hz)) Quantity(T, Snu) {
        if (comptime U_in.dim.eql(Snu_dim)) {
            return input_quantity.to(Snu);
        } else if (comptime U_in.dim.eql(Slamb_dim)) {
            return input_quantity.mulScalar(c).divScalar(freq.pow(2)).to(Snu);
        } else if (comptime U_in.dim.eql(nuSnu_dim)) {
            return input_quantity.divScalar(freq.*).to(Snu);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn toSnuInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), freq: *const Quantity(f64, si.Hz), out: *Quantity(T, Snu)) void {
        if (comptime U_in.dim.eql(Snu_dim)) {
            input_quantity.toInto(Snu, out);
        } else if (comptime U_in.dim.eql(Slamb_dim)) {
            const scalar = c.div(freq.pow(2));
            var temp: Quantity(T, U_in.mul(@TypeOf(scalar).unit)) = .init(out.value);
            input_quantity.mulScalarInto(scalar, &temp);
            temp.toInto(Snu, out);
        } else if (comptime U_in.dim.eql(nuSnu_dim)) {
            var temp: Quantity(T, U_in.div(@TypeOf(freq.*).unit)) = .init(out.value);
            input_quantity.divScalarInto(freq.*, &temp);
            temp.toInto(Snu, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromSnu(comptime T: type, input_quantity: *const Quantity(T, Snu), comptime U_out: Unit, freq: *const Quantity(f64, si.Hz)) Quantity(T, U_out) {
        if (comptime U_out.dim.eql(Snu_dim)) {
            return input_quantity.to(U_out);
        } else if (comptime U_out.dim.eql(Slamb_dim)) {
            return input_quantity.mulScalar(freq.pow(2)).divScalar(c).to(U_out);
        } else if (comptime U_out.dim.eql(nuSnu_dim)) {
            return input_quantity.mulScalar(freq.*).to(U_out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromSnuInto(comptime T: type, input_quantity: *const Quantity(T, Snu), comptime U_out: Unit, freq: *const Quantity(f64, si.Hz), out: *Quantity(T, U_out)) void {
        if (comptime U_out.dim.eql(Snu_dim)) {
            input_quantity.toInto(U_out, out);
        } else if (comptime U_out.dim.eql(Slamb_dim)) {
            const scalar = freq.pow(2).div(c);
            var temp: Quantity(T, Snu.mul(@TypeOf(scalar).unit)) = .init(out.value);
            input_quantity.mulScalarInto(scalar, &temp);
            temp.toInto(U_out, out);
        } else if (comptime U_out.dim.eql(nuSnu_dim)) {
            var temp: Quantity(T, Snu.mul(@TypeOf(freq.*).unit)) = .init(out.value);
            input_quantity.mulScalarInto(freq.*, &temp);
            temp.toInto(U_out, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn toSLnu(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), freq: *const Quantity(f64, si.Hz)) Quantity(T, SLnu) {
        if (comptime U_in.dim.eql(SLnu_dim)) {
            return input_quantity.to(SLnu);
        } else if (comptime U_in.dim.eql(SLlamb_dim)) {
            return input_quantity.mulScalar(c).divScalar(freq.pow(2)).to(SLnu);
        } else if (comptime U_in.dim.eql(nuSLnu_dim)) {
            return input_quantity.divScalar(freq.*).to(SLnu);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn toSLnuInto(comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), freq: *const Quantity(f64, si.Hz), out: *Quantity(T, SLnu)) void {
        if (comptime U_in.dim.eql(SLnu_dim)) {
            input_quantity.toInto(SLnu, out);
        } else if (comptime U_in.dim.eql(SLlamb_dim)) {
            const scalar = c.div(freq.pow(2));
            var temp: Quantity(T, U_in.mul(@TypeOf(scalar).unit)) = .init(out.value);
            input_quantity.mulScalarInto(scalar, &temp);
            temp.toInto(SLnu, out);
        } else if (comptime U_in.dim.eql(nuSLnu_dim)) {
            var temp: Quantity(T, U_in.div(@TypeOf(freq.*).unit)) = .init(out.value);
            input_quantity.divScalarInto(freq.*, &temp);
            temp.toInto(SLnu, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromSLnu(comptime T: type, input_quantity: *const Quantity(T, SLnu), comptime U_out: Unit, freq: *const Quantity(f64, si.Hz)) Quantity(T, U_out) {
        if (comptime U_out.dim.eql(SLnu_dim)) {
            return input_quantity.to(U_out);
        } else if (comptime U_out.dim.eql(SLlamb_dim)) {
            return input_quantity.mulScalar(freq.pow(2)).divScalar(c).to(U_out);
        } else if (comptime U_out.dim.eql(nuSLnu_dim)) {
            return input_quantity.mulScalar(freq.*).to(U_out);
        } else {
            @compileError(unknownConversion(name));
        }
    }

    fn fromSLnuInto(comptime T: type, input_quantity: *const Quantity(T, SLnu), comptime U_out: Unit, freq: *const Quantity(f64, si.Hz), out: *Quantity(T, U_out)) void {
        if (comptime U_out.dim.eql(SLnu_dim)) {
            input_quantity.toInto(U_out, out);
        } else if (comptime U_out.dim.eql(SLlamb_dim)) {
            const scalar = freq.pow(2).div(c);
            var temp: Quantity(T, SLnu.mul(@TypeOf(scalar).unit)) = .init(out.value);
            input_quantity.mulScalarInto(scalar, &temp);
            temp.toInto(U_out, out);
        } else if (comptime U_out.dim.eql(nuSLnu_dim)) {
            var temp: Quantity(T, SLnu.mul(@TypeOf(freq.*).unit)) = .init(out.value);
            input_quantity.mulScalarInto(freq.*, &temp);
            temp.toInto(U_out, out);
        } else {
            @compileError(unknownConversion(name));
        }
    }
};

// Yet to be implemnted:
//
// "beam_angular_area", -> This requires a beam unit defined in astropy in "astrophys". But this is
// a bit misleading in my opinion.
// "logarithmic",
// "thermodynamic_temperature" -> requires cosmology
/// The set of named equivalencies, used as a comptime selector by
/// `Quantity.toWithEquivalency`. Each variant maps to the correspondingly named
/// equivalency struct above.
pub const Equivalency = enum {
    const Self = @This();

    spectral,
    mass_energy,
    temperature_energy,
    temperature,
    molar_mass_amu,
    doppler_redshift,
    dimensionless_angles,
    parallax,
    zero_point_flux,
    doppler_relativistic,
    doppler_optical,
    doppler_radio,
    magnetic_flux_field,
    brightness_temperature,
    plate_scale,
    pixel_scale,
    spectral_density,

    /// Dispatches to the selected equivalency's `convert`. `self` is comptime so the
    /// switch compiles only the chosen branch, which is why equivalencies that take
    /// no `args` never see the (possibly empty) tuple. The `args` tuple is forwarded
    /// only to the equivalencies that require it.
    pub fn convert(comptime self: Self, comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, args: anytype) Quantity(T, U_out) {
        const ArgsType = @TypeOf(args);
        const args_type_info = @typeInfo(ArgsType);
        if (args_type_info != .@"struct") {
            @compileError("Expected tuple or struct argument, found " ++ @typeName(ArgsType));
        }
        return switch (self) {
            .spectral => spectral_equivalency.convert(T, U_in, input_quantity, U_out),
            .mass_energy => mass_energy_equivalency.convert(T, U_in, input_quantity, U_out),
            .temperature_energy => temperature_energy_equivalency.convert(T, U_in, input_quantity, U_out),
            .temperature => temperature_equivalency.convert(T, U_in, input_quantity, U_out),
            .molar_mass_amu => molar_mass_amu_equivalency.convert(T, U_in, input_quantity, U_out),
            .doppler_redshift => doppler_redshift_equivalency.convert(T, U_in, input_quantity, U_out),
            .dimensionless_angles => dimensionless_angles_equivalency.convert(T, U_in, input_quantity, U_out),
            .parallax => parallax_equivalency.convert(T, U_in, input_quantity, U_out),
            .zero_point_flux => zero_point_flux_equivalency.convert(T, U_in, input_quantity, U_out, args),
            .doppler_relativistic => doppler_relativistic_equivalency.convert(T, U_in, input_quantity, U_out, args),
            .doppler_optical => doppler_optical_equivalency.convert(T, U_in, input_quantity, U_out, args),
            .doppler_radio => doppler_radio_equivalency.convert(T, U_in, input_quantity, U_out, args),
            .magnetic_flux_field => magnetic_flux_field_equivalency.convert(T, U_in, input_quantity, U_out, args),
            .brightness_temperature => brightness_temperature_equivalency.convert(T, U_in, input_quantity, U_out, args),
            .plate_scale => plate_scale_equivalency.convert(T, U_in, input_quantity, U_out, args),
            .pixel_scale => pixel_scale_equivalency.convert(T, U_in, input_quantity, U_out, args),
            .spectral_density => spectral_density_equivalency.convert(T, U_in, input_quantity, U_out, args),
        };
    }

    /// Dispatches to the selected equivalency's `convertInto`. `self` is comptime so the
    /// switch compiles only the chosen branch, which is why equivalencies that take
    /// no `args` never see the (possibly empty) tuple. The `args` tuple is forwarded
    /// only to the equivalencies that require it.
    pub fn convertInto(comptime self: Self, comptime T: type, comptime U_in: Unit, input_quantity: *const Quantity(T, U_in), comptime U_out: Unit, out: *Quantity(T, U_out), args: anytype) void {
        const ArgsType = @TypeOf(args);
        const args_type_info = @typeInfo(ArgsType);
        if (args_type_info != .@"struct") {
            @compileError("Expected tuple or struct argument, found " ++ @typeName(ArgsType));
        }
        return switch (self) {
            .spectral => spectral_equivalency.convertInto(T, U_in, input_quantity, U_out, out),
            .mass_energy => mass_energy_equivalency.convertInto(T, U_in, input_quantity, U_out, out),
            .temperature_energy => temperature_energy_equivalency.convertInto(T, U_in, input_quantity, U_out, out),
            .temperature => temperature_equivalency.convertInto(T, U_in, input_quantity, U_out, out),
            .molar_mass_amu => molar_mass_amu_equivalency.convertInto(T, U_in, input_quantity, U_out, out),
            .doppler_redshift => doppler_redshift_equivalency.convertInto(T, U_in, input_quantity, U_out, out),
            .dimensionless_angles => dimensionless_angles_equivalency.convertInto(T, U_in, input_quantity, U_out, out),
            .parallax => parallax_equivalency.convertInto(T, U_in, input_quantity, U_out, out),
            .zero_point_flux => zero_point_flux_equivalency.convertInto(T, U_in, input_quantity, U_out, out, args),
            .doppler_relativistic => doppler_relativistic_equivalency.convertInto(T, U_in, input_quantity, U_out, out, args),
            .doppler_optical => doppler_optical_equivalency.convertInto(T, U_in, input_quantity, U_out, out, args),
            .doppler_radio => doppler_radio_equivalency.convertInto(T, U_in, input_quantity, U_out, out, args),
            .magnetic_flux_field => magnetic_flux_field_equivalency.convertInto(T, U_in, input_quantity, U_out, out, args),
            .brightness_temperature => brightness_temperature_equivalency.convertInto(T, U_in, input_quantity, U_out, out, args),
            .plate_scale => plate_scale_equivalency.convertInto(T, U_in, input_quantity, U_out, out, args),
            .pixel_scale => pixel_scale_equivalency.convertInto(T, U_in, input_quantity, U_out, out, args),
            .spectral_density => spectral_density_equivalency.convertInto(T, U_in, input_quantity, U_out, out, args),
        };
    }
};

test "spectral equivalency" {
    // Wavelength <-> Hz
    const green_wl: Quantity(f64, si.nm) = .init(550);
    const green_freq = spectral_equivalency.convert(f64, si.nm, &green_wl, si.Hz);
    try testing.expectApproxEqRel(5.45077196e14, green_freq.value, 1e-8);
    try testing.expectEqual(si.Hz, green_freq.getUnit());
    const green_wl_back = spectral_equivalency.convert(f64, si.Hz, &green_freq, si.nm);
    try testing.expectApproxEqAbs(green_wl.value, green_wl_back.value, 1e-15);

    // Wavelength <-> Energy
    const lyman: Quantity(f64, si.nm) = .init(121.567);
    const lyman_en = spectral_equivalency.convert(f64, si.nm, &lyman, us.eV);
    try testing.expectApproxEqRel(10.19883673, lyman_en.value, 1e-8);
    try testing.expectEqual(us.eV, lyman_en.getUnit());
    const lyman_back = spectral_equivalency.convert(f64, us.eV, &lyman_en, si.nm);
    try testing.expectApproxEqRel(lyman.value, lyman_back.value, 1e-15);

    // Hz <-> Energy
    const co_line: Quantity(f64, us.GHz) = .init(115.271);
    const co_line_energy = spectral_equivalency.convert(f64, us.GHz, &co_line, si.J);
    try testing.expectApproxEqRel(7.63793732e-23, co_line_energy.value, 1e-8);
    try testing.expectEqual(si.J, co_line_energy.getUnit());
    const co_line_back = spectral_equivalency.convert(f64, si.J, &co_line_energy, us.GHz);
    try testing.expectApproxEqRel(co_line.value, co_line_back.value, 1e-15);

    // Wavelength <-> Wavenumber spectroscopic
    const mir: Quantity(f64, us.um) = .init(10);
    const mir_wn_spec = spectral_equivalency.convert(f64, us.um, &mir, si.cm.pow(-1));
    try testing.expectApproxEqAbs(1000, mir_wn_spec.value, 1e-15);
    try testing.expectEqual(si.cm.pow(-1), mir_wn_spec.getUnit());
    const mir_back = spectral_equivalency.convert(f64, si.cm.pow(-1), &mir_wn_spec, us.um);
    try testing.expectApproxEqRel(mir.value, mir_back.value, 1e-15);

    // Wavelength <-> Wavenumber angular
    const green_wn_ang = spectral_equivalency.convert(f64, us.nm, &green_wl, si.rad.div(si.cm));
    try testing.expectApproxEqRel(114239.73285781, green_wn_ang.value, 1e-8);
    try testing.expectEqual(si.rad.div(si.cm), green_wn_ang.getUnit());
    const green_wl_back_ang = spectral_equivalency.convert(f64, si.rad.div(si.cm), &green_wn_ang, us.nm);
    try testing.expectApproxEqRel(green_wl.value, green_wl_back_ang.value, 1e-15);
}

test "spectral equivalency Into" {
    var buf: [2]f64 = .{ 1, 2 };
    const wl: Quantity([]f64, si.nm) = .init(&buf);
    var buf_res: [2]f64 = undefined;
    var wl_freq: Quantity([]f64, si.Hz) = .init(&buf_res);
    spectral_equivalency.convertInto([]f64, si.nm, &wl, si.Hz, &wl_freq);
    try zatest.expectApproxEqRelIter(@constCast(&[_]f64{ 2.99792458e+17, 1.49896229e+17 }), wl_freq.value, 1e-15);
    try testing.expectEqual(si.Hz, wl_freq.getUnit());
    var buf_back: [2]f64 = undefined;
    var wl_back: Quantity([]f64, si.nm) = .init(&buf_back);
    spectral_equivalency.convertInto([]f64, si.Hz, &wl_freq, si.nm, &wl_back);
    try zatest.expectApproxEqRelIter(wl.value, wl_back.value, 1e-15);
}

test "mass energy equivalency" {
    const m_p_energy = mass_energy_equivalency.convert(f64, si.kg, &cst.m_p.quantity, si.J);
    try testing.expectApproxEqRel(1.50327762e-10, m_p_energy.value, 1e-8);
    try testing.expectEqual(si.J, m_p_energy.getUnit());
    const m_p_back = mass_energy_equivalency.convert(f64, si.J, &m_p_energy, si.kg);
    try testing.expectApproxEqRel(cst.m_p.quantity.value, m_p_back.value, 1e-15);
}

test "mass energy equivalency Into" {
    var buf: [1]f64 = .{10};
    const mass: Quantity([]f64, si.g) = .init(&buf);
    var buf_res: [1]f64 = undefined;
    var mass_energy: Quantity([]f64, us.eV) = .init(&buf_res);
    mass_energy_equivalency.convertInto([]f64, si.g, &mass, us.eV, &mass_energy);
    try zatest.expectApproxEqRelIter(@constCast(&[_]f64{5.6095886e+33}), mass_energy.value, 1e-8);
    try testing.expectEqual(us.eV, mass_energy.getUnit());
    var buf_back: [1]f64 = undefined;
    var mass_back: Quantity([]f64, si.g) = .init(&buf_back);
    mass_energy_equivalency.convertInto([]f64, us.eV, &mass_energy, si.g, &mass_back);
    try zatest.expectApproxEqRelIter(mass.value, mass_back.value, 1e-15);
}

test "temperature energy equivalency" {
    const temp: Quantity(f64, si.K) = .init(300);
    const to_en = temperature_energy_equivalency.convert(f64, si.K, &temp, si.J);
    try testing.expectApproxEqRel(4.141947e-21, to_en.value, 1e-8);
    try testing.expectEqual(si.J, to_en.getUnit());
    const temp_back = temperature_energy_equivalency.convert(f64, si.J, &to_en, si.K);
    try testing.expectApproxEqRel(temp.value, temp_back.value, 1e-15);
}

test "temperature energy equivalency Into" {
    var buf: [2]f64 = .{ 3, 6 };
    const en: Quantity([]f64, si.J) = .init(&buf);
    var buf_res: [2]f64 = undefined;
    var temp: Quantity([]f64, si.K) = .init(&buf_res);
    temperature_energy_equivalency.convertInto([]f64, si.J, &en, si.K, &temp);
    try zatest.expectApproxEqRelIter(@constCast(&[_]f64{ 2.17289115e+23, 4.34578231e+23 }), temp.value, 1e-8);
    try testing.expectEqual(si.K, temp.getUnit());
    var buf_back: [2]f64 = undefined;
    var en_back: Quantity([]f64, si.J) = .init(&buf_back);
    temperature_energy_equivalency.convertInto([]f64, si.K, &temp, si.J, &en_back);
    try zatest.expectApproxEqRelIter(en.value, en_back.value, 1e-15);
}

test "molar mass amu equivalency" {
    const amu: Quantity(f64, misc.u) = .init(1);
    const amu_to_mm = molar_mass_amu_equivalency.convert(f64, misc.u, &amu, si.g.div(si.mol));
    try testing.expectApproxEqRel(1.0, amu_to_mm.value, 1e-8);
    try testing.expectEqual(si.g.div(si.mol), amu_to_mm.getUnit());
    const amu_back = molar_mass_amu_equivalency.convert(f64, si.g.div(si.mol), &amu_to_mm, misc.u);
    try testing.expectApproxEqRel(1.0, amu_back.value, 1e-15);
}

test "molar mass amu equivalency Into" {
    var buf: [1]f64 = .{2};
    const mol_mass: Quantity([]f64, si.kg.div(si.mol)) = .init(&buf);
    var buf_res: [1]f64 = undefined;
    var amu: Quantity([]f64, misc.u) = .init(&buf_res);
    molar_mass_amu_equivalency.convertInto([]f64, @TypeOf(mol_mass).unit, &mol_mass, misc.u, &amu);
    try zatest.expectApproxEqRelIter(@constCast(&[_]f64{2000}), amu.value, 1e-8);
    try testing.expectEqual(misc.u, amu.getUnit());
    var buf_back: [1]f64 = undefined;
    var mol_mass_back: Quantity([]f64, si.kg.div(si.mol)) = .init(&buf_back);
    molar_mass_amu_equivalency.convertInto([]f64, misc.u, &amu, si.kg.div(si.mol), &mol_mass_back);
    try zatest.expectApproxEqAbsIter(mol_mass.value, mol_mass_back.value, 1e-15);
}

test "doppler redshift equivalency" {
    const cmb_z: Quantity(f64, UNITLESS) = .init(1089);
    const cmb_v = doppler_redshift_equivalency.convert(f64, UNITLESS, &cmb_z, us.km.div(si.s));
    try testing.expectApproxEqRel(299791.9533418, cmb_v.value, 1e-8);
    try testing.expectEqual(us.km.div(si.s), cmb_v.getUnit());
    const cmb_z_back = doppler_redshift_equivalency.convert(f64, us.km.div(si.s), &cmb_v, UNITLESS);
    try testing.expectApproxEqRel(cmb_z.value, cmb_z_back.value, 1e-8);
}

test "doppler redshift equivalency Into" {
    var buf: [2]f64 = .{ 1, 2 };
    const z: Quantity([]f64, UNITLESS) = .init(&buf);
    var buf_res: [2]f64 = undefined;
    var v: Quantity([]f64, us.km.div(si.s)) = .init(&buf_res);
    doppler_redshift_equivalency.convertInto([]f64, UNITLESS, &z, us.km.div(si.s), &v);
    try zatest.expectApproxEqRelIter(@constCast(&[_]f64{ 179875.4748, 239833.9664 }), v.value, 1e-8);
    try testing.expectEqual(us.km.div(si.s), v.getUnit());
    var buf_back: [2]f64 = undefined;
    var z_back: Quantity([]f64, UNITLESS) = .init(&buf_back);
    doppler_redshift_equivalency.convertInto([]f64, us.km.div(si.s), &v, UNITLESS, &z_back);
    try zatest.expectApproxEqRelIter(z.value, z_back.value, 1e-15);
}

test "dimensionless angles equivalency" {
    const pi_rad: Quantity(f64, si.rad) = .init(std.math.pi);
    const pi_dimless = dimensionless_angles_equivalency.convert(f64, si.rad, &pi_rad, UNITLESS);
    try testing.expectApproxEqAbs(std.math.pi, pi_dimless.value, 1e-15);
    try testing.expectEqual(UNITLESS, pi_dimless.getUnit());
    const pi_rad_back = dimensionless_angles_equivalency.convert(f64, UNITLESS, &pi_dimless, si.rad);
    try testing.expectApproxEqAbs(std.math.pi, pi_rad_back.value, 1e-15);
    try testing.expectEqual(si.rad, pi_rad_back.getUnit());

    const one_rad: Quantity(f64, si.rad) = .init(1);
    const one_dimless = dimensionless_angles_equivalency.convert(f64, si.rad, &one_rad, UNITLESS);
    try testing.expectApproxEqAbs(1.0, one_dimless.value, 1e-15);
    const one_rad_back = dimensionless_angles_equivalency.convert(f64, UNITLESS, &one_dimless, si.rad);
    try testing.expectApproxEqAbs(1.0, one_rad_back.value, 1e-15);

    const deg_180: Quantity(f64, si.deg) = .init(180);
    const deg_180_dimless = dimensionless_angles_equivalency.convert(f64, si.deg, &deg_180, UNITLESS);
    try testing.expectApproxEqAbs(std.math.pi, deg_180_dimless.value, 1e-15);
    try testing.expectEqual(UNITLESS, deg_180_dimless.getUnit());
    const deg_180_back = dimensionless_angles_equivalency.convert(f64, UNITLESS, &deg_180_dimless, si.deg);
    try testing.expectApproxEqAbs(180, deg_180_back.value, 1e-15);
    try testing.expectEqual(si.deg, deg_180_back.getUnit());
}

test "dimensionless angles equivalency Into" {
    var buf: [2]f64 = .{ 1, std.math.pi };
    const rad: Quantity([]f64, si.rad) = .init(&buf);
    var buf_res: [2]f64 = undefined;
    var rad_dimless: Quantity([]f64, UNITLESS) = .init(&buf_res);
    dimensionless_angles_equivalency.convertInto([]f64, si.rad, &rad, UNITLESS, &rad_dimless);
    try zatest.expectApproxEqAbsIter(@constCast(&[_]f64{ 1, std.math.pi }), rad_dimless.value, 1e-15);
    try testing.expectEqual(UNITLESS, rad_dimless.getUnit());
    var buf_back: [2]f64 = undefined;
    var rad_back: Quantity([]f64, si.rad) = .init(&buf_back);
    dimensionless_angles_equivalency.convertInto([]f64, UNITLESS, &rad_dimless, si.rad, &rad_back);
    try zatest.expectApproxEqAbsIter(buf, rad_back.value, 1e-15);
}

test "parallax equivalency" {
    const as_1: Quantity(f64, si.arcsec) = .init(1);
    const pc = parallax_equivalency.convert(f64, si.arcsec, &as_1, us.pc);
    try testing.expectApproxEqAbs(1.0, pc.value, 1e-15);
    try testing.expectEqual(us.pc, pc.getUnit());
    const as_1_back = parallax_equivalency.convert(f64, us.pc, &pc, si.arcsec);
    try testing.expectApproxEqAbs(1.0, as_1_back.value, 1e-15);
}

test "parallax equivalency Into" {
    var buf: [2]f64 = .{ 2, 30 };
    const pc: Quantity([]f64, us.pc) = .init(&buf);
    var buf_res: [2]f64 = undefined;
    var as: Quantity([]f64, us.arcsec) = .init(&buf_res);
    parallax_equivalency.convertInto([]f64, us.pc, &pc, us.arcsec, &as);
    try zatest.expectApproxEqAbsIter(@constCast(&[_]f64{ 0.5, 0.033333333333333 }), as.value, 1e-15);
    try testing.expectEqual(us.arcsec, as.getUnit());
    var buf_back: [2]f64 = undefined;
    var pc_back: Quantity([]f64, us.pc) = .init(&buf_back);
    parallax_equivalency.convertInto([]f64, us.arcsec, &as, us.pc, &pc_back);
    try zatest.expectApproxEqAbsIter(buf, pc_back.value, 1e-15);
}

test "zero point flux equivalency" {
    const mag: Quantity(f64, photo.mgy) = .init(1);
    const flux0: Quantity(f64, us.Jy) = .init(3631);
    const mag_to_flux = zero_point_flux_equivalency.convert(f64, photo.mgy, &mag, us.Jy, .{flux0});
    try testing.expectApproxEqAbs(3631, mag_to_flux.value, 1e-15);
    try testing.expectEqual(us.Jy, mag_to_flux.getUnit());
    const mag_back = zero_point_flux_equivalency.convert(f64, us.Jy, &mag_to_flux, photo.mgy, .{flux0});
    try testing.expectApproxEqAbs(1, mag_back.value, 1e-15);
    try testing.expectEqual(photo.mgy, mag_back.getUnit());
}

test "zero point flux equivalency Into" {
    var buf: [2]f64 = .{ 2, 5 };
    const f: Quantity([]f64, us.Jy) = .init(&buf);
    const f0: Quantity(f64, us.Jy) = .init(100);
    var buf_res: [2]f64 = undefined;
    var buf_back: [2]f64 = undefined;
    var res: Quantity([]f64, us.mgy) = .init(&buf_res);
    zero_point_flux_equivalency.convertInto([]f64, us.Jy, &f, us.mgy, &res, .{f0});
    try zatest.expectApproxEqAbsIter(@constCast(&[_]f64{ 0.02, 0.05 }), res.value, 1e-15);
    try testing.expectEqual(us.mgy, res.getUnit());
    var f_back: Quantity([]f64, us.Jy) = .init(&buf_back);
    zero_point_flux_equivalency.convertInto([]f64, us.mgy, &res, us.Jy, &f_back, .{f0});
    try zatest.expectApproxEqAbsIter(f.value, f_back.value, 1e-15);
}

test "doppler relativistic equivalency" {
    const v_unit: Unit = comptime us.km.div(si.s);
    const f0: Quantity(f64, us.MHz) = .init(1420.405751768);
    const v1000: Quantity(f64, v_unit) = .init(1000);
    const v1000_f = doppler_relativistic_equivalency.convert(f64, v_unit, &v1000, us.MHz, .{f0});
    try testing.expectApproxEqAbs(1415.67566395, v1000_f.value, 1e-8);
    try testing.expectEqual(us.MHz, v1000_f.getUnit());
    const v1000_back = doppler_relativistic_equivalency.convert(f64, us.MHz, &v1000_f, v_unit, .{f0});
    try testing.expectApproxEqRel(1000, v1000_back.value, 1e-8);

    const v_09c: Quantity(f64, v_unit) = .init(270_000);
    const v_09c_f = doppler_relativistic_equivalency.convert(f64, v_unit, &v_09c, us.MHz, .{f0});
    try testing.expectApproxEqAbs(324.79343469, v_09c_f.value, 1e-8);
    const v_09c_back = doppler_relativistic_equivalency.convert(f64, us.MHz, &v_09c_f, v_unit, .{f0});
    try testing.expectApproxEqRel(270_000, v_09c_back.value, 1e-15);

    const wl_21: Quantity(f64, si.cm) = .init(21.106);
    const wl_21f_v = doppler_relativistic_equivalency.convert(f64, si.cm, &wl_21, v_unit, .{f0});
    try testing.expectApproxEqAbs(-1.62003622, wl_21f_v.value, 1e-8);
}

test "doppler relativistic equivalency Into" {
    const f0: Quantity(f64, us.GHz) = .init(90);

    var buf: [2]f64 = .{ 45, 30 };
    const freq: Quantity([]f64, us.GHz) = .init(&buf);
    var buf_res: [2]f64 = undefined;
    var v: Quantity([]f64, us.km.div(us.s)) = .init(&buf_res);
    doppler_relativistic_equivalency.convertInto([]f64, us.GHz, &freq, us.km.div(us.s), &v, .{f0});
    try zatest.expectApproxEqRelIter(@constCast(&[_]f64{ 179875.4748, 239833.9664 }), v.value, 1e-15);
    try testing.expectEqual(v.getUnit(), us.km.div(us.s));
    var buf_back: [2]f64 = undefined;
    var freq_back: Quantity([]f64, us.GHz) = .init(&buf_back);
    doppler_relativistic_equivalency.convertInto([]f64, us.km.div(us.s), &v, us.GHz, &freq_back, .{f0});
    try zatest.expectApproxEqRelIter(buf_back, buf, 1e-15);
}

test "doppler optical equivalency" {
    const v_unit: Unit = comptime us.km.div(si.s);
    const f0: Quantity(f64, us.MHz) = .init(1420.405751768);
    const v_1000: Quantity(f64, v_unit) = .init(1000);
    const v_1000_f = doppler_optical_equivalency.convert(f64, v_unit, &v_1000, us.MHz, .{f0});
    try testing.expectApproxEqAbs(1415.6835397776715, v_1000_f.value, 1e-8);
    try testing.expectEqual(us.MHz, v_1000_f.getUnit());
    const v_1000_back = doppler_optical_equivalency.convert(f64, us.MHz, &v_1000_f, v_unit, .{f0});
    try testing.expectApproxEqRel(1000, v_1000_back.value, 1e-8);
}

test "doppler optical equivalency Into" {
    const f0: Quantity(f64, us.GHz) = .init(90);

    const freq: Quantity([]f64, us.GHz) = .init(@constCast(&[_]f64{ 56.25, 50 }));
    var v_buf: [2]f64 = undefined;
    var v: Quantity([]f64, us.km.div(si.s)) = .init(&v_buf);
    doppler_optical_equivalency.convertInto([]f64, us.GHz, &freq, us.km.div(us.s), &v, .{f0});
    try zatest.expectApproxEqRelIter(@constCast(&[_]f64{ 179875.4748, 239833.9664 }), v.value, 1e-15);
    try testing.expectEqual(us.km.div(us.s), v.getUnit());
    var buf_back: [2]f64 = undefined;
    var freq_back: Quantity([]f64, us.GHz) = .init(&buf_back);
    doppler_optical_equivalency.convertInto([]f64, us.km.div(us.s), &v, us.GHz, &freq_back, .{f0});
    try zatest.expectApproxEqRelIter(buf_back, freq.value, 1e-15);
}

test "doppler radio equivalency" {
    const v_unit: Unit = comptime us.km.div(si.s);
    const f0: Quantity(f64, us.MHz) = .init(1420.405751768);
    const v_1000: Quantity(f64, v_unit) = .init(1000);
    const v_1000_f = doppler_radio_equivalency.convert(f64, v_unit, &v_1000, us.MHz, .{f0});
    try testing.expectApproxEqAbs(1415.66778817, v_1000_f.value, 1e-8);
    try testing.expectEqual(us.MHz, v_1000_f.getUnit());
    const v_1000_back = doppler_radio_equivalency.convert(f64, us.MHz, &v_1000_f, v_unit, .{f0});
    try testing.expectApproxEqRel(1000, v_1000_back.value, 1e-8);
}

test "doppler radio equivalency Into" {
    const f0: Quantity(f64, us.GHz) = .init(90);

    const freq: Quantity([]f64, us.GHz) = .init(@constCast(&[_]f64{ 36, 18 }));
    var v_buf: [2]f64 = undefined;
    var v: Quantity([]f64, us.km.div(si.s)) = .init(&v_buf);
    doppler_radio_equivalency.convertInto([]f64, us.GHz, &freq, us.km.div(us.s), &v, .{f0});
    try zatest.expectApproxEqRelIter(@constCast(&[_]f64{ 179875.4748, 239833.9664 }), v.value, 1e-15);
    try testing.expectEqual(us.km.div(us.s), v.getUnit());
    var buf_back: [2]f64 = undefined;
    var freq_back: Quantity([]f64, us.GHz) = .init(&buf_back);
    doppler_radio_equivalency.convertInto([]f64, us.km.div(us.s), &v, us.GHz, &freq_back, .{f0});
    try zatest.expectApproxEqRelIter(buf_back, freq.value, 1e-15);
}

test "magnetic flux field equivalency" {
    const mu_r: Quantity(f64, UNITLESS) = .init(1);
    const H: Quantity(f64, us.Oe) = .init(1);
    const B_G = magnetic_flux_field_equivalency.convert(f64, us.Oe, &H, us.G, .{mu_r});
    try testing.expectApproxEqAbs(1.0, B_G.value, 1e-8);
    const H_from_B_G = magnetic_flux_field_equivalency.convert(f64, us.G, &B_G, us.Oe, .{mu_r});
    try testing.expectApproxEqRel(H.value, H_from_B_G.value, 1e-15);

    const B_T = magnetic_flux_field_equivalency.convert(f64, us.Oe, &H, us.T, .{mu_r});
    try testing.expectApproxEqAbs(1e-4, B_T.value, 1e-8);
    const H_from_B_T = magnetic_flux_field_equivalency.convert(f64, us.T, &B_T, us.Oe, .{mu_r});
    try testing.expectApproxEqRel(H.value, H_from_B_T.value, 1e-15);
}

test "magnetic flux field equivalency Into" {
    const H: Quantity([]f64, us.Oe) = .init(@constCast(&[_]f64{12}));
    var buf_res: [1]f64 = undefined;
    var B_T: Quantity([]f64, us.T) = .init(&buf_res);
    magnetic_flux_field_equivalency.convertInto([]f64, us.Oe, &H, us.T, &B_T, .{});
    try zatest.expectApproxEqRelIter(@constCast(&[_]f64{12e-4}), B_T.value, 1e-8);

    const mu_r: Quantity(f64, UNITLESS) = .init(24);
    var buf_back: [1]f64 = undefined;
    var H_back: Quantity([]f64, us.Oe) = .init(&buf_back);
    magnetic_flux_field_equivalency.convertInto([]f64, us.T, &B_T, us.Oe, &H_back, .{mu_r});
    try zatest.expectApproxEqRelIter(@constCast(&[_]f64{0.5}), H_back.value, 1e-15);
}

test "brightness temperature equivalency" {
    // TODO: Check properly the reason why increasing the branch quota is needed.
    // This could be improve or solved, maybe.
    @setEvalBranchQuota(2000);
    const wb: Quantity(f64, si.arcsec.pow(2)) = .init(std.math.pi * 2 * 2500);
    const jy_q: Quantity(f64, us.Jy) = .init(1);
    const jy_sr_q = jy_q.div(wb).to(us.Jy.div(si.sr));
    const nu: Quantity(f64, us.GHz) = .init(5);
    const as_K = brightness_temperature_equivalency.convert(f64, us.Jy.div(si.sr), &jy_sr_q, si.K, .{nu});
    try testing.expectApproxEqAbs(3.52629392, as_K.value, 1e-8);
    const jy_sr_back = brightness_temperature_equivalency.convert(f64, si.K, &as_K, us.Jy.div(si.sr), .{nu});
    try testing.expectApproxEqRel(jy_sr_q.value, jy_sr_back.value, 1e-15);

    // VLA
    const bmin: Quantity(f64, si.arcsec) = .init(15);
    const bmaj = bmin;
    const fwhm_to_sig: Quantity(f64, UNITLESS) = .init(1.0 / @sqrt(8.0 * @log(2.0)));
    const wb_2: Quantity(f64, si.sr) = fwhm_to_sig.pow(2).mul(bmin).mul(bmaj).mulValue(2 * std.math.pi).to(si.sr);
    const jans: Quantity(f64, us.Jy) = .init(1);
    const as_K_2 = brightness_temperature_equivalency.convert(f64, us.Jy, &jans, si.K, .{ nu, wb_2 });
    try testing.expectApproxEqAbs(217.26583888, as_K_2.value, 1e-8);
    const jans_back = brightness_temperature_equivalency.convert(f64, si.K, &as_K_2, us.Jy, .{ nu, wb_2 });
    try testing.expectApproxEqRel(jans.value, jans_back.value, 1e-15);
}

test "brightness temperature equivalency Into" {
    const nu: Quantity(f64, us.GHz) = .init(1.4);

    const sb: Quantity([]f64, us.Jy.div(si.sr)) = .init(@constCast(&[_]f64{ 1.0, 2.0 }));
    var buf_res: [2]f64 = undefined;
    var res: Quantity([]f64, si.K) = .init(&buf_res);
    brightness_temperature_equivalency.convertInto([]f64, us.Jy.div(si.sr), &sb, si.K, &res, .{nu});
    try zatest.expectApproxEqRelIter(@constCast(&[_]f64{ 1.66062685e-05, 3.32125370e-05 }), res.value, 1e-8);

    var buf_back: [2]f64 = undefined;
    var back: Quantity([]f64, us.Jy.div(si.sr)) = .init(&buf_back);
    brightness_temperature_equivalency.convertInto([]f64, si.K, &res, us.Jy.div(si.sr), &back, .{nu});
    try zatest.expectApproxEqRelIter(sb.value, back.value, 1e-8);
}

test "plate scale equivalency" {
    const angle: Quantity(f64, si.arcsec) = .init(1);
    const plate_scale: Quantity(f64, si.m.div(si.rad)) = .init(10);
    const as_length = plate_scale_equivalency.convert(f64, si.arcsec, &angle, si.m, .{plate_scale});
    try testing.expectApproxEqAbs(4.84813681e-05, as_length.value, 1e-8);
    const angle_back = plate_scale_equivalency.convert(f64, si.m, &as_length, si.arcsec, .{plate_scale});
    try testing.expectApproxEqRel(angle.value, angle_back.value, 1e-15);
}

test "plate scale equivalency Into" {
    const angle: Quantity([]f64, si.deg) = .init(@constCast(&[_]f64{ 10, 15 }));
    const plate_scale: Quantity(f64, si.m.div(si.rad)) = .init(1);
    var buf_res: [2]f64 = undefined;
    var length: Quantity([]f64, si.m) = .init(&buf_res);
    plate_scale_equivalency.convertInto([]f64, si.deg, &angle, si.m, &length, .{plate_scale});
    try zatest.expectApproxEqAbsIter(@constCast(&[2]f64{ 0.174532925199, 0.261799387799 }), length.value, 1e-12);
    var buf_back: [2]f64 = undefined;
    var angle_back: Quantity([]f64, si.deg) = .init(&buf_back);
    plate_scale_equivalency.convertInto([]f64, si.m, &length, si.deg, &angle_back, .{plate_scale});
    try zatest.expectApproxEqRelIter(angle.value, angle_back.value, 1e-15);
}

test "pixel scale equivalency" {
    {
        const pixscale: Quantity(f64, si.arcsec.div(misc.pix)) = .init(0.4);
        const pix: Quantity(f64, misc.pix) = .init(240);
        const as_angle = pixel_scale_equivalency.convert(f64, misc.pix, &pix, si.arcmin, .{pixscale});
        try testing.expectApproxEqAbs(1.6, as_angle.value, 1e-8);
        const pix_back = pixel_scale_equivalency.convert(f64, si.arcmin, &as_angle, misc.pix, .{pixscale});
        try testing.expectApproxEqRel(240, pix_back.value, 1e-15);
    }
    {
        const pixscale: Quantity(f64, misc.pix.div(us.cm)) = .init(100);
        const pix: Quantity(f64, misc.pix) = .init(1024);
        const as_cm = pixel_scale_equivalency.convert(f64, misc.pix, &pix, si.cm, .{pixscale});
        try testing.expectApproxEqAbs(10.24, as_cm.value, 1e-8);
        const pix_back = pixel_scale_equivalency.convert(f64, si.cm, &as_cm, misc.pix, .{pixscale});
        try testing.expectApproxEqRel(pix.value, pix_back.value, 1e-15);
    }
}

test "pixel scale equivalency Into" {
    const pixscale: Quantity(f64, si.arcsec.div(misc.pix)) = .init(0.2);

    const px: Quantity([]f64, misc.pix) = .init(@constCast(&[_]f64{ 10, 25 }));
    var ang_buf: [2]f64 = undefined;
    var ang: Quantity([]f64, us.arcsec) = .init(&ang_buf);
    pixel_scale_equivalency.convertInto([]f64, misc.pix, &px, us.arcsec, &ang, .{pixscale});
    try zatest.expectApproxEqRelIter(@constCast(&[_]f64{ 2, 5 }), ang.value, 1e-12);

    var px_buf: [2]f64 = undefined;
    var px_back: Quantity([]f64, misc.pix) = .init(&px_buf);
    pixel_scale_equivalency.convertInto([]f64, us.arcsec, &ang, misc.pix, &px_back, .{pixscale});
    try zatest.expectApproxEqRelIter(@constCast(&[_]f64{ 10, 25 }), px_back.value, 1e-12);
}

test "spectral density equivalency" {
    const nu_w: Quantity(f64, si.AA) = .init(3500);
    const spe = spectral_density_equivalency;
    // Flux
    {
        const f: Quantity(f64, us.Jy) = .init(1.5);
        const as_nuFnu = spe.convert(f64, us.Jy, &f, si.W.div(si.m.pow(2)), .{nu_w});
        try testing.expectApproxEqRel(1.28482482e-11, as_nuFnu.value, 1e-8);
        const f_back = spe.convert(f64, si.W.div(si.m.pow(2)), &as_nuFnu, us.Jy, .{nu_w});
        try testing.expectApproxEqRel(f.value, f_back.value, 1e-15);
    }
    // Luminosity density
    {
        const L: Quantity(f64, si.W.div(si.Hz)) = .init(1.5);
        const as_pow = spe.convert(f64, si.W.div(si.Hz), &L, si.W, .{nu_w});
        try testing.expectApproxEqRel(1.28482482e+15, as_pow.value, 1e-8);
        const L_back = spe.convert(f64, si.W, &as_pow, si.W.div(si.Hz), .{nu_w});
        try testing.expectApproxEqRel(L.value, L_back.value, 1e-15);
    }
    // Surface brightness
    {
        const sb: Quantity(f64, spe.Snu) = .init(1.5);
        const as_nuSnu = spe.convert(f64, spe.Snu, &sb, spe.nuSnu, .{nu_w});
        try testing.expectApproxEqRel(1.28482482e+15, as_nuSnu.value, 1e-8);
        const sb_back = spe.convert(f64, spe.nuSnu, &as_nuSnu, spe.Snu, .{nu_w});
        try testing.expectApproxEqRel(sb.value, sb_back.value, 1e-15);
    }
    // Luminosity surface brightness
    {
        const lsb: Quantity(f64, spe.SLnu) = .init(1.5);
        const as_SLlamb = spe.convert(f64, spe.SLnu, &lsb, spe.SLlamb, .{nu_w});
        try testing.expectApproxEqRel(3.67092806e+21, as_SLlamb.value, 1e-8);
        const lsb_back = spe.convert(f64, spe.SLlamb, &as_SLlamb, spe.SLnu, .{nu_w});
        try testing.expectApproxEqRel(lsb.value, lsb_back.value, 1e-15);
    }
}

test "spectral density equivalency Into" {
    const nu_w: Quantity(f64, si.AA) = .init(3500);
    const W_m2 = comptime si.W.div(si.m.pow(2));
    const L_Hz = comptime si.W.div(si.Hz);
    const spe = spectral_density_equivalency;

    // Flux
    {
        const f: Quantity([]f64, us.Jy) = .init(@constCast(&[_]f64{ 1.5, 3.0 }));
        var out_buf: [2]f64 = undefined;
        var as_nuFnu: Quantity([]f64, W_m2) = .init(&out_buf);
        spe.convertInto([]f64, us.Jy, &f, W_m2, &as_nuFnu, .{nu_w});
        try zatest.expectApproxEqRelIter(@constCast(&[_]f64{ 1.28482482e-11, 2.56964964e-11 }), as_nuFnu.value, 1e-8);

        var buf_back: [2]f64 = undefined;
        var f_back: Quantity([]f64, us.Jy) = .init(&buf_back);
        spe.convertInto([]f64, W_m2, &as_nuFnu, us.Jy, &f_back, .{nu_w});
    }

    // Luminosity density
    {
        const L: Quantity([]f64, L_Hz) = .init(@constCast(&[_]f64{ 1.5, 3.0 }));
        var out_buf: [2]f64 = undefined;
        var as_pow: Quantity([]f64, si.W) = .init(&out_buf);
        spe.convertInto([]f64, L_Hz, &L, si.W, &as_pow, .{nu_w});
        try zatest.expectApproxEqRelIter(@constCast(&[_]f64{ 1.28482482e+15, 2.56964964e+15 }), as_pow.value, 1e-8);

        var back_buf: [2]f64 = undefined;
        var L_back: Quantity([]f64, L_Hz) = .init(&back_buf);
        spe.convertInto([]f64, si.W, &as_pow, L_Hz, &L_back, .{nu_w});
        try zatest.expectApproxEqRelIter(@constCast(&[_]f64{ 1.5, 3.0 }), L_back.value, 1e-12);
    }

    // Surface brightness
    {
        const sb: Quantity([]f64, spe.Snu) = .init(@constCast(&[_]f64{ 1.5, 3.0 }));
        var out_buf: [2]f64 = undefined;
        var as_nuSnu: Quantity([]f64, spe.nuSnu) = .init(&out_buf);
        spe.convertInto([]f64, spe.Snu, &sb, spe.nuSnu, &as_nuSnu, .{nu_w});
        try zatest.expectApproxEqRelIter(@constCast(&[_]f64{ 1.28482482e+15, 2.56964964e+15 }), as_nuSnu.value, 1e-8);

        var back_buf: [2]f64 = undefined;
        var sb_back: Quantity([]f64, spe.Snu) = .init(&back_buf);
        spe.convertInto([]f64, spe.nuSnu, &as_nuSnu, spe.Snu, &sb_back, .{nu_w});
        try zatest.expectApproxEqRelIter(@constCast(&[_]f64{ 1.5, 3.0 }), sb_back.value, 1e-12);
    }

    // Luminosity surface brightness
    {
        const lsb: Quantity([]f64, spe.SLnu) = .init(@constCast(&[_]f64{ 1.5, 3.0 }));
        var out_buf: [2]f64 = undefined;
        var as_SLlamb: Quantity([]f64, spe.SLlamb) = .init(&out_buf);
        spe.convertInto([]f64, spe.SLnu, &lsb, spe.SLlamb, &as_SLlamb, .{nu_w});
        try zatest.expectApproxEqRelIter(@constCast(&[_]f64{ 3.67092806e+21, 7.34185612e+21 }), as_SLlamb.value, 1e-8);

        var back_buf: [2]f64 = undefined;
        var lsb_back: Quantity([]f64, spe.SLnu) = .init(&back_buf);
        spe.convertInto([]f64, spe.SLlamb, &as_SLlamb, spe.SLnu, &lsb_back, .{nu_w});
        try zatest.expectApproxEqRelIter(@constCast(&[_]f64{ 1.5, 3.0 }), lsb_back.value, 1e-12);
    }
}
