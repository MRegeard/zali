const std = @import("std");
const unit = @import("unit.zig");
const Unit = unit.Unit;
const utils = @import("utils.zig");
const si = @import("si.zig");
const cst = @import("../constants.zig");
const defUnitFromUnit = utils.defUnitFromUnit;
const defUnitFromConst = utils.defUnitFromConst;

const PI_F64: f64 = std.math.pi;


// Areas
//
pub const barn: Unit = defUnitFromUnit(si.m.pow(2), 1e-28, .initUniqueSymbol("barn"));

// Angular measurements
//
pub const cycle: Unit = defUnitFromUnit(si.sr, 2.0 * PI_F64, .initUniqueSymbol("cycle"));
pub const spat: Unit = defUnitFromUnit(si.sr, 4.0 * PI_F64, .initUniqueSymbol("spat"));

// Pressure
//
pub const bar: Unit = defUnitFromUnit(si.Pa, 1e5, .initUniqueSymbol("bar"));
// Astropy comment:
// The Torr is almost the same as mmHg but not quite.
// See https://en.wikipedia.org/wiki/Torr
pub const Torr: Unit = defUnitFromUnit(si.Pa, cst.atm.quantity.value / 760, .initUniqueSymbol("Torr"));

// Mass
//
pub const M_p: Unit = defUnitFromConst(cst.m_p, .initUniqueSymbol("M_p"));
pub const M_e: Unit = defUnitFromConst(cst.m_e, .initUniqueSymbol("M_e"));
pub const u: Unit = defUnitFromConst(cst.u, .initUniqueSymbol("u"));

// Energy
//
pub const eV: Unit = defUnitFromUnit(si.J, cst.e.quantity.value, .initUniqueSymbol("eV"));
pub const Ry: Unit = defUnitFromUnit(eV, cst.Ryd.quantity.to(eV).value, .initUniqueSymbol("Ry"));

// Computer
//
pub const bit: Unit = defUnitFromUnit(unit.UNITLESS, 1.0, .initUniqueSymbol("bit"));
pub const byte: Unit = defUnitFromUnit(bit, 8, .initUniqueSymbol("byte"));
pub const pix: Unit = defUnitFromUnit(unit.UNITLESS, 1.0, .initUniqueSymbol("pix"));
pub const vox: Unit = defUnitFromUnit(unit.UNITLESS, 1.0, .initUniqueSymbol("vox"));
