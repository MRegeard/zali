const std = @import("std");
const cst = @import("../constants.zig");
const utils = @import("utils.zig");
const unit = @import("unit.zig");
const si = @import("si.zig");
const Unit = unit.Unit;
const defUnitFromUnit = utils.defUnitFromUnit;
const defUnitFromConst = utils.defUnitFromConst;

const PI_F64: f64 = std.math.pi;

// Length
//
pub const au: Unit = defUnitFromConst(cst.au, .initUniqueSymbol("au"));
pub const pc: Unit = defUnitFromConst(cst.pc, .initUniqueSymbol("pc"));
pub const kpc: Unit = defUnitFromConst(cst.kpc, .initUniqueSymbol("kpc"));
pub const solRad: Unit = defUnitFromConst(cst.R_sun, .initUniqueSymbol("solRad"));
pub const jupiterRad: Unit = defUnitFromConst(cst.R_jup, .initUniqueSymbol("jupiterRad"));
pub const earthRad: Unit = defUnitFromConst(cst.R_earth, .initUniqueSymbol("earthRad"));
pub const lyr: Unit = defUnitFromUnit(cst.c.quantity.mul(si.yr).to(si.m), 1.0, .initUniqueSymbol("lyr"));
pub const lsec: Unit = defUnitFromUnit(cst.c.mul(si.s).to(si.m), 1.0, .initUniqueSymbol("lsec"));

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
