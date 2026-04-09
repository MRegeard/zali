const std = @import("std");
const testing = std.testing;
const unit = @import("unit.zig");
const Unit = unit.Unit;
const UNITLESS = unit.UNITLESS;
const quantityMod = @import("quantity.zig");
const Quantity = quantityMod.Quantity;
const isQuantity = quantityMod.isQuantity;
const si = @import("si.zig");
const photo = @import("photometric.zig");
const misc = @import("misc.zig");
const cst = @import("../constants.zig");
const imperial = @import("imperial.zig");

pub const ConversionError = error{
    UnknownSpectralConversion,
    UnknownZeroPointFluxConversion,
    UnknownMassEnergyConversion,
    UnknownTemperatureEnergyConversion,
    UnknownMolarMassAmuConversion,
    UnknownDopplerRedshiftConversion,
    UnknownDopplerRelativisticConversion,
    UnknownDopplerOpticalConversion,
    UnknownTemperatureConversion,
    UnknownMagneticFluxFieldConversion,
    UnknownDopplerRadioConversion,
    UnknownPixelScaleConversion,
    UnexpectedPixelPower,
    UnexpectedConversionArgument,
    ExpectedConversionArgument,
};

const InOutUnit = struct {
    const Self = @This();

    in: Unit,
    out: Unit,

    pub fn cmp(self: *const Self, other_in: Unit, other_out: Unit) bool {
        return (self.in.dim.eql(other_in.dim) and self.out.dim.eql(other_out.dim));
    }

    pub fn eqlExact(self: *const Self, other_in: Unit, other_out: Unit) bool {
        return (self.in.eqlExact(other_in) and self.out.eqlExact(other_out));
    }
};

pub fn spectralConvert(comptime T: type, comptime U_in: Unit, input_quantity: Quantity(T, U_in), U_out: Unit) ConversionError!Quantity(T, U_out) {
    const two_pi: Quantity(f64, unit.UNITLESS) = .init(std.math.pi);
    const inv_m_spec: Unit = si.m.pow(-1);
    const inv_m_ang: Unit = si.rad.div(si.m);
    const iou = InOutUnit{ .in = U_in, .out = U_out };
    // m    -> HZ
    if (iou.cmp(si.m, si.Hz)) {
        const to_m = try input_quantity.to(si.m);
        const to_hz = try to_m.pow(-1).mul(cst.c.quantity).to(si.Hz);
        return try to_hz.to(U_out);
    }
    // m -> J
    if (iou.cmp(si.m, si.J)) {
        const to_m = try input_quantity.to(si.m);
        const to_J = try to_m.pow(-1).mul(cst.c.quantity).mul(cst.h.quantity).to(si.J);
        return try to_J.to(U_out);
    }
    // Hz -> J
    if (iou.cmp(si.Hz, si.J)) {
        const to_hz = try input_quantity.to(si.Hz);
        const to_J = try to_hz.mul(cst.h.quantity).to(si.J);
        return try to_J.to(U_out);
    }
    // J -> Hz
    if (iou.cmp(si.J, si.Hz)) {
        const to_J = try input_quantity.to(si.J);
        const to_hz = try to_J.div(cst.h.quantity).to(si.Hz);
        return try to_hz.to(U_out);
    }
    // m -> m-1(spec)
    if (iou.cmp(si.m, inv_m_spec)) {
        const to_m_inv = input_quantity.pow(-1).to(inv_m_spec);
        return try to_m_inv.to(U_out);
    }
    // Hz -> m-1(spec)
    if (iou.cmp(si.Hz, inv_m_spec)) {
        const to_Hz = try input_quantity.to(si.Hz);
        const to_m_inv = try to_Hz.div(cst.c.quantity).to(inv_m_spec);
        return try to_m_inv.to(U_out);
    }
    // m-1(spec) -> Hz
    if (iou.cmp(inv_m_spec, si.Hz)) {
        const to_m_inv = try input_quantity.to(inv_m_spec);
        const to_Hz = try to_m_inv.mul(cst.c.quantity).to(si.Hz);
        return try to_Hz.to(U_out);
    }
    // J -> m-1(spec)
    if (iou.cmp(si.J, inv_m_spec)) {
        const to_J = try input_quantity.to(si.J);
        const to_m_inv = try to_J.div(cst.c.quantity.mul(cst.h.quantity)).to(inv_m_spec);
        return try to_m_inv.to(U_out);
    }
    // m-1(spec) -> J
    if (iou.cmp(inv_m_spec, si.J)) {
        const to_m_inv = try input_quantity.to(inv_m_spec);
        const to_J = try to_m_inv.mul(cst.c.quantity.mul(cst.h.quantity)).to(si.J);
        return try to_J.to(U_out);
    }
    // m-1(spec) -> m-1(ang)
    if (iou.cmp(inv_m_spec, inv_m_ang)) {
        const to_m_inv_spec = try input_quantity.to(inv_m_spec);
        const to_m_inv_ang = try to_m_inv_spec.mul(two_pi).to(inv_m_ang);
        return try to_m_inv_ang.to(U_out);
    }
    // m-1(ang) -> m-1(spec)
    if (iou.cmp(inv_m_ang, inv_m_spec)) {
        const to_m_inv_ang = try input_quantity.to(inv_m_ang);
        const to_m_inv_spec = try to_m_inv_ang.div(two_pi).to(inv_m_spec);
        return try to_m_inv_spec.to(U_out);
    }
    // m -> m-1(ang)
    if (iou.cmp(si.m, inv_m_ang)) {
        const to_m = try input_quantity.to(si.m);
        const to_m_inv_ang = try two_pi.div(to_m).to(inv_m_ang);
        return try to_m_inv_ang.to(U_out);
    }
    // Hz -> m-1(ang)
    if (iou.cmp(si.Hz, inv_m_ang)) {
        const to_Hz = try input_quantity.to(si.Hz);
        const to_m_inv_ang = try to_Hz.mul(two_pi).div(cst.c.quantity).to(inv_m_ang);
        return try to_m_inv_ang.to(U_out);
    }
    // m-1(ang) -> Hz
    if (iou.cmp(inv_m_ang, si.Hz)) {
        const to_m_inv_ang = try input_quantity.to(inv_m_ang);
        const to_Hz = try to_m_inv_ang.mul(cst.c.quantity).div(two_pi);
        return try to_Hz.to(U_out);
    }
    // J -> m-1(ang)
    if (iou.cmp(si.J, inv_m_ang)) {
        const to_J = try input_quantity.to(si.J);
        const to_m_inv_ang = try to_J.mul(two_pi).div(cst.h.quantity.mul(cst.c.quantity)).to(inv_m_ang);
        return try to_m_inv_ang.to(U_out);
    }
    // m-1(ang) -> J
    if (iou.cmp(inv_m_ang, si.J)) {
        const to_m_inv_ang = try input_quantity.to(inv_m_ang);
        const to_J = try to_m_inv_ang.mul(cst.h.quantity.mul(cst.c.quantity)).div(two_pi).to(si.J);
        return try to_J.to(U_out);
    }
    return ConversionError.UnknownSpectralConversion;
}

pub fn zeroPointFluxConvert(comptime T: type, comptime U_in: Unit, input_quantity: Quantity(T, U_in), comptime U_out: Unit, args: anytype) ConversionError!Quantity(T, U_out) {
    if (!U_in.eqlExact(photo.mgy)) {
        @compileError("Expected Quantity in units of maggy for ZeroPointFlux equivalency");
    }
    const args_fields = @typeInfo(@TypeOf(args)).@"struct".fields;
    if (args_fields.len != 1) {
        @compileError("Expected one arguments for ZeroPointFlux equivalency");
    }
    const f_type = args_fields[0].type;
    const f_type_info = @typeInfo(f_type);
    // TODO: check that it is a Quantity type. For now it is assumed.
    if (f_type_info == .@"struct" and f_type.unit.eqlExact(U_out)) {
        return try input_quantity.mul(@field(args, args_fields[0].name)).to(U_out);
    }
    return ConversionError.UnknownZeroPointFluxConversion;
}

pub fn massEnergyConvert(comptime T: type, comptime U_in: Unit, input_quantity: Quantity(T, U_in), comptime U_out: Unit) ConversionError!Quantity(T, U_out) {
    const c2 = cst.c.quantity.pow(2);
    const c2_splat: Quantity(T, c2.getUnit()) = .initScalarValue(c2);
    const iou = InOutUnit{ .in = U_in, .out = U_out };
    if (iou.cmp(si.kg, si.J)) {
        const to_kg = try input_quantity.to(si.kg);
        const to_J = try to_kg.mul(c2_splat).to(si.J);
        return try to_J.to(U_out);
    }
    if (iou.cmp(si.J, si.kg)) {
        const to_J = try input_quantity.to(si.J);
        const to_kg = try to_J.div(c2_splat).to(si.kg);
        return try to_kg.to(U_out);
    }
    if (iou.cmp(si.kg.div(si.m.pow(2)), si.J.div(si.m.pow(2)))) {
        const to_input_unit = try input_quantity.to(si.kg.div(si.m.pow(2)));
        const to_output_unit = try to_input_unit.mul(c2_splat).to(si.J.div(si.m.pow(2)));
        return try to_output_unit.to(U_out);
    }
    if (iou.cmp(si.J.div(si.m.pow(2)), si.kg.div(si.m.pow(2)))) {
        const to_input_unit = try input_quantity.to(si.J.div(si.m.pow(2)));
        const to_output_unit = try to_input_unit.mul(c2_splat).to(si.kg.div(si.m.pow(2)));
        return try to_output_unit.to(U_out);
    }
    if (iou.cmp(si.kg.div(si.m.pow(3)), si.J.div(si.m.pow(3)))) {
        const to_input_unit = try input_quantity.to(si.kg.div(si.m.pow(3)));
        const to_output_unit = try to_input_unit.mul(c2_splat).to(si.J.div(si.m.pow(3)));
        return try to_output_unit.to(U_out);
    }
    if (iou.cmp(si.J.div(si.m.pow(3)), si.kg.div(si.m.pow(3)))) {
        const to_input_unit = try input_quantity.to(si.J.div(si.m.pow(3)));
        const to_output_unit = try to_input_unit.mul(c2_splat).to(si.kg.div(si.m.pow(3)));
        return try to_output_unit.to(U_out);
    }
    if (iou.cmp(si.kg.div(si.s), si.J.div(si.s))) {
        const to_input_unit = try input_quantity.to(si.kg.div(si.s));
        const to_output_unit = try to_input_unit.mul(c2_splat).to(si.J.div(si.s));
        return try to_output_unit.to(U_out);
    }
    if (iou.cmp(si.J.div(si.s), si.kg.div(si.s))) {
        const to_input_unit = try input_quantity.to(si.J.div(si.s));
        const to_output_unit = try to_input_unit.div(c2_splat).to(si.kg.div(si.s));
        return try to_output_unit.to(U_out);
    }
    return ConversionError.UnknownMassEnergyConversion;
}

pub fn temperatureEnergyConvert(comptime T: type, comptime U_in: Unit, input_quantity: Quantity(T, U_in), U_out: Unit) ConversionError!Quantity(T, U_out) {
    const iou = InOutUnit{ .in = U_in, .out = U_out };
    if (iou.cmp(si.K, misc.eV)) {
        const to_K = try input_quantity.to(si.K);
        const to_eV = try to_K.mul(cst.k_B.quantity).to(misc.eV);
        return try to_eV.to(U_out);
    }
    if (iou.cmp(misc.eV, si.K)) {
        const to_eV = try input_quantity.to(misc.eV);
        const to_K = try to_eV.div(cst.k_B.quantity).to(si.K);
        return try to_K.to(U_out);
    }
    return ConversionError.UnknownTemperatureEnergyConversion;
}

pub fn molarMassAmuConvert(comptime T: type, comptime U_in: Unit, input_quantity: Quantity(T, U_in), U_out: Unit) ConversionError!Quantity(T, U_out) {
    const iou = InOutUnit{ .in = U_in, .out = U_out };
    if (iou.cmp(si.g.div(si.mol), misc.u)) {
        const to_g_per_mol = try input_quantity.to(si.g.div(si.mol));
        const to_u: Quantity(T, misc.u) = .init(to_g_per_mol.value);
        return try to_u.to(U_out);
    }
    if (iou.cmp(misc.u, si.g.div(si.mol))) {
        const to_u = try input_quantity.to(misc.u);
        const to_g_per_mol: Quantity(T, si.g.div(si.mol)) = .init(to_u.value);
        return try to_g_per_mol.to(U_out);
    }
    return ConversionError.UnknownMolarMassAmuConversion;
}

pub fn dopplerRedshiftConvert(comptime T: type, comptime U_in: Unit, input_quantity: Quantity(T, U_in), U_out: Unit) ConversionError!Quantity(T, U_out) {
    const iou = InOutUnit{ .in = U_in, .out = U_out };
    const rv_unit: Unit = si.km.div(si.s);
    const c_kms = cst.c.quantity.to(rv_unit);
    const one: Quantity(T, UNITLESS) = .initScalarValue(1);
    if (iou.cmp(UNITLESS, rv_unit)) {
        const zponesq = input_quantity.add(one).pow(2);
        return try c_kms.mul(zponesq.sub(one)).div(zponesq.add(one)).to(U_out);
    }
    if (iou.cmp(rv_unit, UNITLESS)) {
        const rv: Quantity(T, rv_unit) = .initScalarValue(1);
        const beta = rv.div(c_kms);
        return try beta.add(one).div(one.sub(beta)).sqrt().sub(one).to(U_out);
    }
    return ConversionError.UnknownDopplerRedshiftConversion;
}

// TODO check that args is a quantity
// TODO check the whole function because it is messy ! (and late right now !)
pub fn dopplerRelativisticConvert(comptime T: type, comptime U_in: Unit, input_quantity: Quantity(T, U_in), U_out: Unit, args: anytype) ConversionError!Quantity(T, U_out) {
    const args_fields = @typeInfo(@TypeOf(args)).@"struct".fields;
    if (args_fields.len != 1) {
        @compileError("Expected one arguments for ZeroPointFlux equivalency");
    }
    // TODO: check that it is a Quantity type. For now it is assumed.
    const rest = @field(args, args_fields[0].name);
    assertIsSpectralunit(rest);
    const iou = InOutUnit{ .in = U_in, .out = U_out };
    const rv_unit: Unit = si.km.div(si.s);
    const c_kms = cst.c.quantity.to(rv_unit);
    const one: Quantity(T, UNITLESS) = .initScalarValue(1);
    if (iou.cmp(si.Hz, rv_unit)) {
        const restfreq = try rest.toWithEquivalency(si.Hz, .Spectral, .{});
        const restfreq2 = restfreq.pow(2);
        const iq_Hz = try input_quantity.to(si.Hz);
        const iq2 = iq_Hz.pow(2);
        return try restfreq2.sub(iq2).div(restfreq2.add(iq2)).mul(c_kms).to(U_out);
    }
    if (iou.cmp(rv_unit, si.Hz)) {
        const iq_kms = try input_quantity.to(rv_unit);
        const voverc = iq_kms.div(c_kms);
        const restfreq = try rest.toWithEquivalency(si.Hz, .Spectral, .{});
        const voverc_p1 = one.add(voverc);
        const voverc_s1 = one.sub(voverc);
        return try restfreq.mul(voverc_s1.div(voverc_p1).sqrt()).to(U_out);
    }
    if (iou.cmp(si.AA, rv_unit)) {
        const restwav = try rest.toWithEquivalency(si.AA, .Spectral, .{});
        const restwav2 = restwav.pow(2);
        const iq_AA = try input_quantity.to(si.AA);
        const iq2 = iq_AA.pow(2);
        return try iq2.sub(restwav2).div(restwav2.add(iq2)).mul(c_kms);
    }
    if (iou.cmp(rv_unit, si.AA)) {
        const iq_kms = try input_quantity.to(rv_unit);
        const voverc = iq_kms.div(c_kms);
        const voverc_p1 = one.add(voverc);
        const voverc_s1 = one.sub(voverc);
        const restwav = try rest.toWithEquivalency(si.AA, .Spectral, .{});
        return try restwav.mul(voverc_p1.div(voverc_s1).sqrt()).to(U_out);
    }
    if (iou.cmp(misc.eV, rv_unit)) {
        const resten = try rest.toWithEquivalency(misc.eV, .Spectral, .{});
        const resten2 = resten.pow(2);
        const iq_eV = try input_quantity.to(misc.eV);
        const iq2 = iq_eV.pow(2);
        return try resten2.sub(iq2).div(resten2.add(iq2)).mul(c_kms).to(U_out);
    }
    if (iou.cmp(rv_unit, misc.eV)) {
        const iq_kms = try input_quantity.to(rv_unit);
        const voverc = iq_kms.div(c_kms);
        const voverc_p1 = one.add(voverc);
        const voverc_s1 = one.sub(voverc);
        const resten = try rest.toWithEquivalency(misc.eV, .Spectral, .{});
        return try resten.mul(voverc_s1.div(voverc_p1).sqrt()).to(U_out);
    }
    return ConversionError.UnknownDopplerRelativisticConversion;
}

pub fn dopplerOpticalConvert(comptime T: type, comptime U_in: Unit, input_quantity: Quantity(T, U_in), U_out: Unit, args: anytype) ConversionError!Quantity(T, U_out) {
    const args_fields = @typeInfo(@TypeOf(args)).@"struct".fields;
    if (args_fields.len != 1) {
        @compileError("Expected one arguments for ZeroPointFlux equivalency");
    }
    // TODO: check that it is a Quantity type. For now it is assumed.
    const rest = @field(args, args_fields[0].name);
    assertIsSpectralunit(rest);
    const iou = InOutUnit{ .in = U_in, .out = U_out };
    const rv_unit: Unit = si.km.div(si.s);
    const c_kms = cst.c.quantity.to(rv_unit);
    const one: Quantity(T, UNITLESS) = .initScalarValue(1);
    if (iou.cmp(misc.eV, rv_unit)) {
        const iq_eV = try input_quantity.to(misc.eV);
        const resten = try rest.toWithEquivalency(misc.eV, .Spectral, .{});
        const to_rv = resten.div(iq_eV).sub(one).mul(c_kms);
        return try to_rv.to(U_out);
    }
    if (iou.cmp(rv_unit, misc.eV)) {
        const iq_rv = try input_quantity.to(rv_unit);
        const iq_p1 = iq_rv.add(one);
        const resten = try rest.toWithEquivalency(misc.eV, .Spectral, .{});
        const to_eV = resten.div(iq_p1.mul(c_kms));
        return try to_eV.to(U_out);
    }
    if (iou.cmp(si.Hz, rv_unit)) {
        const iq_Hz = try input_quantity.to(si.Hz);
        const restfreq = try rest.toWithEquivalency(si.Hz, .Spectral, .{});
        const to_rv = restfreq.div(iq_Hz).sub(one).mul(c_kms);
        return try to_rv.to(U_out);
    }
    if (iou.cmp(rv_unit, si.Hz)) {
        const iq_rv = try input_quantity.to(rv_unit);
        const iq_p1 = iq_rv.add(one);
        const restfreq = try rest.toWithEquivalency(si.Hz, .Spectral, .{});
        const to_Hz = restfreq.div(iq_p1.div(c_kms));
        return try to_Hz.to(U_out);
    }
    if (iou.cmp(si.AA, rv_unit)) {
        const iq_AA = try input_quantity.to(si.AA);
        const restwav = try rest.toWithEquivalency(si.AA, .Spectral, .{});
        const to_wav = iq_AA.div(restwav).sub(one).mul(c_kms);
        return try to_wav.to(U_out);
    }
    if (iou.cmp(rv_unit, si.AA)) {
        const iq_rv = try input_quantity.to(rv_unit);
        const iq_p1 = iq_rv.add(one);
        const restwav = try rest.toWithEquivalency(si.AA, .Spectral, .{});
        const to_rv = restwav.mul(iq_p1.div(c_kms));
        return try to_rv.to(U_out);
    }
    return ConversionError.UnknownDopplerOpticalConverion;
}

pub fn temperatureConvert(comptime T: type, comptime U_in: Unit, input_quantity: Quantity(T, U_in), comptime U_out: Unit) ConversionError!Quantity(T, U_out) {
    const iou = InOutUnit{ .in = U_in, .out = U_out };
    const factor: Quantity(T, UNITLESS) = .initScalarValue(1.8);
    // Note: For each conversion, the given offset is written in the destination unit such that
    // we can add to it.
    // We create a new quantity in the destination unit initialised with the input_quantity value
    // and then we add/sub/mul/div the relevant quantity to conevert from the input unit to the
    // destination unit.
    if (iou.eqlExact(si.K, si.degC)) {
        // C = K - 273.15
        const offset_K_degC: Quantity(T, si.degC) = .initScalarValue(273.15);
        const in_C: Quantity(T, si.degC) = .init(input_quantity.value);
        return in_C.sub(offset_K_degC);
    }
    if (iou.eqlExact(si.degC, si.K)) {
        // K = C + 273.15
        const offset_K_degC: Quantity(T, si.K) = .initScalarValue(273.15);
        const in_K: Quantity(T, si.K) = .init(input_quantity.value);
        return in_K.add(offset_K_degC);
    }
    if (iou.eqlExact(si.degC, imperial.degF)) {
        // F = C * 1.8 + 32
        const offset_degC_degF: Quantity(T, imperial.degF) = .initScalarValue(32);
        const in_degF: Quantity(T, imperial.degF) = .init(input_quantity.value);
        return in_degF.mul(factor).add(offset_degC_degF);
    }
    if (iou.eqlExact(imperial.degF, si.degC)) {
        // C = (F - 32) / 1.8
        const offset_degC_degF: Quantity(T, si.degC) = .initScalarValue(32);
        const in_degC: Quantity(T, si.degC) = .init(input_quantity.value);
        return in_degC.sub(offset_degC_degF).div(factor);
    }
    if (iou.eqlExact(si.K, imperial.degF)) {
        // F = K * 1.8 - 459.67
        const offset_K_degF: Quantity(T, imperial.degF) = .initScalarValue(459.67);
        const in_degF: Quantity(T, imperial.degF) = .init(input_quantity.value);
        return in_degF.mul(factor).sub(offset_K_degF);
    }
    if (iou.eqlExact(imperial.degF, si.K)) {
        // K = (F + 459.67) / 1.8
        const offset_K_degF: Quantity(T, si.K) = .initScalarValue(459.67);
        const in_K: Quantity(T, si.K) = .init(input_quantity.value);
        return in_K.add(offset_K_degF).div(factor);
    }
    return ConversionError.UnknownTemperatureConversion;
}

pub fn magneticFluxFieldConvert(comptime T: type, comptime U_in: Unit, input_quantity: Quantity(T, U_in), comptime U_out: Unit, args: anytype) ConversionError!Quantity(T, U_out) {
    var mu_r: Quantity(T, UNITLESS) = .initScalarValue(1);
    const args_fields = @typeInfo(@TypeOf(args)).@"struct".fields;
    if (args_fields.len == 1) {
        mu_r.value = @field(args, args_fields[0].name).value;
    } else if (args_fields.len > 1) {
        @compileError("Expected one arguments for ZeroPointFlux equivalency");
    }
    // TODO: check that it is a Quantity type. For now it is assumed.
    // TODO: Add the possibility to input scalar quantity as mu_r is a dimensionless quantity ?
    const iou = InOutUnit{ .in = U_in, .out = U_out };
    if (iou.cmp(si.T, si.A.div(si.m))) {
        const to_T = try input_quantity.to(si.T);
        const to_dest_unit = try to_T.div(cst.mu0.quantity.mul(mu_r)).to(si.A.div(si.m));
        return try to_dest_unit.to(U_out);
    }
    if (iou.cmp(si.A.div(si.m), si.T)) {
        const to_in_unit = try input_quantity.to(si.A.div(si.m));
        const to_T = try to_in_unit.mul(cst.mu0.quantity).mul(mu_r);
        return try to_T.to(U_out);
    }
    return ConversionError.UnknownMagneticFluxFieldConversion;
}

pub fn dopplerRadioConvert(comptime T: type, comptime U_in: Unit, input_quantity: Quantity(T, U_in), U_out: Unit, args: anytype) ConversionError!Quantity(T, U_out) {
    const args_fields = @typeInfo(@TypeOf(args)).@"struct".fields;
    if (args_fields.len != 1) {
        @compileError("Expected one arguments for ZeroPointFlux equivalency");
    }
    // TODO: check that it is a Quantity type. For now it is assumed.
    const rest = @field(args, args_fields[0].name);
    assertIsSpectralunit(rest);
    const iou = InOutUnit{ .in = U_in, .out = U_out };
    const rv_unit: Unit = si.km.div(si.s);
    const c_kms = cst.c.quantity.to(rv_unit);
    const one: Quantity(T, UNITLESS) = .initScalarValue(1);
    if (iou.cmp(misc.eV, rv_unit)) {
        const resten = try rest.toWithEquivalency(misc.eV, .Spectral, .{});
        const iq_eV = try input_quantity.to(misc.eV);
        const to_rv = one.sub(iq_eV.div(resten)).mul(c_kms);
        return try to_rv.to(U_out);
    }
    if (iou.cmp(rv_unit, misc.eV)) {
        const resten = try rest.toWithEquivalency(misc.eV, .Spectral, .{});
        const iq_rv = try input_quantity.to(rv_unit);
        const to_eV = resten.mul(one.sub(iq_rv.div(c_kms)));
        return try to_eV.to(U_out);
    }
    if (iou.cmp(si.Hz, rv_unit)) {
        const restfreq = try rest.toWithEquivalency(si.Hz, .Spectral, .{});
        const iq_Hz = try input_quantity.to(si.Hz);
        const to_rv = one.sub(iq_Hz.div(restfreq)).mul(c_kms);
        return try to_rv.to(U_out);
    }
    if (iou.cmp(rv_unit, si.Hz)) {
        const restfreq = try rest.toWithEquivalency(si.Hz, .Spectral, .{});
        const iq_rv = try input_quantity.to(rv_unit);
        const to_Hz = restfreq.mul(one.sub(iq_rv.div(c_kms)));
        return try to_Hz.to(U_out);
    }
    if (iou.cmp(si.AA, rv_unit)) {
        const restwav = try rest.toWithEquivalency(si.AA, .Spectral, .{});
        const iq_AA = try input_quantity.to(si.AA);
        const to_rv = one.sub(restwav.div(iq_AA)).mul(c_kms);
        return try to_rv.to(U_out);
    }
    if (iou.cpm(rv_unit, si.AA)) {
        const restwav = try rest.toWithEquivalency(si.AA, .Spectral, .{});
        const iq_rv = try input_quantity.to(rv_unit);
        const to_AA = restwav.div(one.sub(iq_rv.div(c_kms)));
        return try to_AA.to(U_out);
    }
    return ConversionError.UnknownDopplerRadioConversion;
}

// TODO: I am a bit confused here, check that !
pub fn pixelScaleConvert(comptime T: type, comptime U_in: Unit, input_quantity: Quantity(T, U_in), U_out: Unit, args: anytype) ConversionError!Quantity(T, U_out) {
    const args_fields = @typeInfo(@TypeOf(args)).@"struct".fields;
    if (args_fields.len != 1) {
        @compileError("Expected one arguments for ZeroPointFlux equivalency");
    }
    // TODO: check that it is a Quantity type. For now it is assumed.
    const pixscale = @field(args, args_fields[0].name);
    const pixscale_symbol = pixscale.getUnit().symbol;
    if (pixscale_symbol.indexOfSymbol(misc.pix.symbol.terms[0])) |idx| {
        const pix_power = pixscale_symbol.terms[idx].exponent;
        var physical_unit: Unit = undefined;
        if (pix_power.eqlScalar(-1)) {
            physical_unit = pixscale.getUnit().mul(misc.pix);
        } else if (pix_power.eqlScalar(1)) {
            physical_unit = pixscale.getUnit().div(misc.pix);
        } else return ConversionError.UnexpectedPixelPower;
        if (input_quantity.getUnit().eqlExact(misc.pix) and physical_unit.eql(U_out)) {
            const pixscale_unitless: Quantity(T, UNITLESS) = .init(pixscale.value);
            const q_out: Quantity(T, physical_unit) = .init(input_quantity.value);
            return try q_out.mul(pixscale_unitless).to(U_out);
        }
        if (input_quantity.getUnit().eql(physical_unit) and misc.pix.eqlExact(U_out)) {
            const pixscale_unitless: Quantity(T, UNITLESS) = .init(pixscale.value);
            const q_out: Quantity(T, misc.pix) = .init(input_quantity.value);
            return q_out.mul(pixscale_unitless);
        }
    }
    return ConversionError.UnknownPixelScaleConversion;
}

// Yet to be implemnted:
//
// RAAD: Require Angle As Dimension
//
// "beam_angular_area", -> RAAD
// "brightness_temperature", -> RAAD
// "dimensionless_angles", -> RAAD
// "logarithmic",
// "parallax", -> RAAD
// "plate_scale",
// "spectral_density",
// "thermodynamic_temperature" -> requires cosmology
pub const Equivalency = enum {
    const Self = @This();

    Spectral,
    ZeroPointFlux,
    MassEnergy,
    TemperatureEnergy,
    MolaMassAmu,
    DopplerRedshift,
    DopplerRelativistic,
    DopplerOptical,
    Temperature,
    MagneticFluxField,
    DopplerRadio,
    PixelScale,

    pub fn convert(self: Self, comptime T: type, comptime U_in: Unit, input_quantity: Quantity(T, U_in), comptime U_out: Unit, args: anytype) ConversionError!Quantity(T, U_out) {
        const ArgsType = @TypeOf(args);
        const args_type_info = @typeInfo(ArgsType);
        if (args_type_info != .@"struct") {
            @compileError("Expected tuple or struct argument, found " ++ @typeName(ArgsType));
        }
        switch (self) {
            .Spectral => {
                if (isEmptyTuple(&args)) {
                    return try spectralConvert(T, U_in, input_quantity, U_out);
                }
                return ConversionError.UnexpectedConversionArgument;
            },
            .ZeroPointFLux => {
                if (!isEmptyTuple(&args)) {
                    return try zeroPointFluxConvert(T, U_in, input_quantity, U_out);
                }
                return ConversionError.ExpectedConversionArgument;
            },
            .MassEnergy => {
                if (isEmptyTuple(&args)) {
                    return try massEnergyConvert(T, U_in, input_quantity, U_out);
                }
                return ConversionError.UnexpectedConversionArgument;
            },
            .TemperatureEnergy => {
                if (isEmptyTuple(&args)) {
                    return try temperatureEnergyConvert(T, U_in, input_quantity, U_out);
                }
                return ConversionError.UnexpectedConversionArgument;
            },
            .MolaMassAmu => {
                if (isEmptyTuple(&args)) {
                    return try molarMassAmuConvert(T, U_in, input_quantity, U_out);
                }
                return ConversionError.UnexpectedConversionArgument;
            },
            .DopplerRedshift => {
                if (isEmptyTuple(&args)) {
                    return try dopplerRedshiftConvert(T, U_in, input_quantity, U_out);
                }
                return ConversionError.UnexpectedConversionArgument;
            },
            .DopplerRelativistic => {
                if (!isEmptyTuple(&args)) {
                    return try dopplerRelativisticConvert(T, U_in, input_quantity, U_out, args);
                }
                return ConversionError.ExpectedConversionArgument;
            },
            .DopplerOptical => {
                if (!isEmptyTuple(&args)) {
                    return try dopplerOpticalConvert(T, U_in, input_quantity, U_out, args);
                }
                return ConversionError.ExpectedConversionArgument;
            },
            .Temperature => {
                if (isEmptyTuple(&args)) {
                    return try temperatureConvert(T, U_in, input_quantity, U_out);
                }
                return ConversionError.UnexpectedConversionArgument;
            },
            .MagneticFluxField => {
                // Note: Here we can accept either no args since mu_r defaults to 1.0
                // or mu_r as an args.
                return try magneticFluxFieldConvert(T, U_in, input_quantity, U_out, args);
            },
            .DopplerRadio => {
                if (!isEmptyTuple(&args)) {
                    return try dopplerRadioConvert(T, U_in, input_quantity, U_out, args);
                }
                return ConversionError.ExpectedConversionArgument;
            },
            .PixelScale => {
                if (!isEmptyTuple(&args)) {
                    return try pixelScaleConvert(T, U_in, input_quantity, U_out, args);
                }
                return ConversionError.ExpectedConversionArgument;
            },
        }
    }
};

fn isEmptyTuple(args: anytype) bool {
    const args_type_info = @typeInfo(@TypeOf(args));
    switch (args_type_info) {
        .pointer => |ptr| return isEmptyTuple(ptr.*),
        .@"struct" => |s| {
            if (s.fields.len == 0 and s.decls.len == 0 and s.is_tuple == true) return true;
            return false;
        },
        else => return false,
    }
}

fn assertIsSpectralunit(quantity: anytype) void {
    if (isQuantity(quantity)) {
        quantity.toWithEquivalency(si.Hz, .Spectral) catch @panic("`args` must be a spectral equivalent (frequency, wavelength or energy).");
    }
}

test "isEmptyTuple" {
    try testing.expect(isEmptyTuple(.{}) == true);
    try testing.expect(isEmptyTuple(&.{}) == true);
    try testing.expect(isEmptyTuple(.{ 12, 45 }) == false);
}
