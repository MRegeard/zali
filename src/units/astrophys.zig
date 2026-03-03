const std = @import("std");
const cst = @import("../constants.zig");
const utils = @import("utils.zig");
const unit = @import("unit.zig");
const si = @import("si.zig");
const Unit = unit.Unit;
const Quantity = @import("quantity.zig").Quantity;
const defUnitFromUnit = utils.defUnitFromUnit;
const defUnitFromConst = utils.defUnitFromConst;

const PI_F64: f64 = std.math.pi;
const yr_quantity: Quantity(f64, si.yr) = .init(1.0);
const s_quantity: Quantity(f64, si.s) = .init(1.0);
const lyr_quantity = cst.c.quantity.mul(yr_quantity).to(si.m) catch unreachable;
const lsec_quantity = cst.c.quantity.mul(s_quantity).to(si.m) catch unreachable;

// Length
//
pub const au: Unit = defUnitFromConst(cst.au, .initUniqueSymbol("au"));
pub const pc: Unit = defUnitFromConst(cst.pc, .initUniqueSymbol("pc"));
pub const kpc: Unit = defUnitFromConst(cst.kpc, .initUniqueSymbol("kpc"));
pub const solRad: Unit = defUnitFromConst(cst.R_sun, .initUniqueSymbol("solRad"));
pub const jupiterRad: Unit = defUnitFromConst(cst.R_jup, .initUniqueSymbol("jupiterRad"));
pub const earthRad: Unit = defUnitFromConst(cst.R_earth, .initUniqueSymbol("earthRad"));
pub const lyr: Unit = defUnitFromUnit(si.m, lyr_quantity.value, .initUniqueSymbol("lyr"));
pub const lsec: Unit = defUnitFromUnit(si.m, lsec_quantity.value, .initUniqueSymbol("lsec"));

// Mass
//
pub const solMass: Unit = defUnitFromConst(cst.M_sun, .initUniqueSymbol("solMass"));
pub const jupiterMass: Unit = defUnitFromConst(cst.M_jup, .initUniqueSymbol("jupiterMass"));
pub const earthMass: Unit = defUnitFromConst(cst.M_earth, .initUniqueSymbol("earthMass"));

// Energy
//
pub const foe: Unit = defUnitFromUnit(si.g.mul(si.cm.pow(2)).div(si.s.pow(2)), 1e51, .initUniqueSymbol("foe"));

// Illumination
//
pub const solLum: Unit = defUnitFromConst(cst.L_sun, .initUniqueSymbol("solLum"));

// Spectral density
//
pub const ph: Unit = defUnitFromUnit(unit.UNITLESS, 1.0, .initUniqueSymbol("ph"));
pub const Jy: Unit = defUnitFromUnit(si.W.div(si.m.pow(2)).div(si.Hz), 1e-26, .initUniqueSymbol("Jy"));
pub const R: Unit = defUnitFromUnit(ph.mul(si.m.pow(-2)).mul(si.s.pow(-1)).mul(si.sr.pow(-1)), 1e10 / (4 * PI_F64), .initUniqueSymbol("R"));

// Events
//
pub const ct: Unit = defUnitFromUnit(unit.UNITLESS, 1.0, .initUniqueSymbol("ct"));
pub const adu: Unit = defUnitFromUnit(unit.UNITLESS, 1.0, .initUniqueSymbol("adu"));
pub const DN: Unit = defUnitFromUnit(unit.UNITLESS, 1.0, .initUniqueSymbol("DN"));
