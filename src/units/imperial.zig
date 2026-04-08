/// Imperial units
///
/// Those are not re-exported to zali.units, and are therefore only available via
/// zali.units.imperial
const std = @import("st");
const si = @import("si.zig");
const unit = @import("unit.zig");
const utils = @import("utils.zig");
const dimMod = @import("dim.zig");
const Unit = unit.Unit;
const defUnitFromUnit = utils.defUnitFromUnit;

// Length
pub const inch: Unit = defUnitFromUnit(si.cm, 2.54, .initUniqueSymbol("inch"));
pub const ft: Unit = defUnitFromUnit(inch, 12, .initUniqueSymbol("ft"));
pub const yd: Unit = defUnitFromUnit(ft, 3, .initUniqueSymbol("yd"));
pub const mi: Unit = defUnitFromUnit(ft, 5280, .initUniqueSymbol("mi"));
pub const mil: Unit = defUnitFromUnit(inch, 0.001, .initUniqueSymbol("mil"));
pub const nmi: Unit = defUnitFromUnit(si.m, 1852, .initUniqueSymbol("nmi"));
pub const fur: Unit = defUnitFromUnit(ft, 660, .initUniqueSymbol("fur"));

// Area
pub const ac: Unit = defUnitFromUnit(ft.pow(2), 43560, .initUniqueSymbol("ac"));

// Volume
pub const gallon: Unit = defUnitFromUnit(si.l, 1 / 0.264172052, .initUniqueSymbol("gallon"));
pub const quart: Unit = defUnitFromUnit(gallon, 1 / 4, .initUniqueSymbol("quart"));
pub const pint: Unit = defUnitFromUnit(quart, 1 / 2, .initUniqueSymbol("pint"));
pub const cup: Unit = defUnitFromUnit(pint, 1 / 2, .initUniqueSymbol("cup"));
pub const foz: Unit = defUnitFromUnit(cup, 1 / 8, .initUniqueSymbol("foz"));
pub const tbsp: Unit = defUnitFromUnit(foz, 1 / 2, .initUniqueSymbol("tbsp"));
pub const tsp: Unit = defUnitFromUnit(tbsp, 1 / 3, .initUniqueSymbol("tsp"));

// Mass
pub const oz: Unit = defUnitFromUnit(si.g, 28.349523125, .initUniqueSymbol("oz"));
pub const lb: Unit = defUnitFromUnit(oz, 16, .initUniqueSymbol("lb"));
pub const st: Unit = defUnitFromUnit(lb, 14, .initUniqueSymbol("st"));
pub const ton: Unit = defUnitFromUnit(lb, 2000, .initUniqueSymbol("ton"));
pub const slug: Unit = defUnitFromUnit(lb, 32.174049, .initUniqueSymbol("slug"));

// Speed
pub const kn: Unit = defUnitFromUnit(nmi.div(si.h), 1, .initUniqueSymbol("kn"));

// Force
pub const lbf: Unit = defUnitFromUnit(slug.mul(ft).mul(si.s.pow(-2)), 1, .initUniqueSymbol("lbf"));
pub const kip: Unit = defUnitFromUnit(lbf, 1000, .initUniqueSymbol("kip"));

// Energy
pub const BTU: Unit = defUnitFromUnit(si.J, 1055.05585, .initUniqueSymbol("BTU"));
pub const cal: Unit = defUnitFromUnit(si.J, 4.184, .initUniqueSymbol("cal"));
pub const kcal: Unit = defUnitFromUnit(cal, 1000, .initUniqueSymbol("kcal"));

// Pressure
pub const psi: Unit = defUnitFromUnit(lbf.mul(inch.pow(-2)), 1, .initUniqueSymbol("psi"));

// Power
pub const hp: Unit = defUnitFromUnit(si.W, 0.00134102209, .initUniqueSymbol("hp"));

// Temperature
pub const degF: Unit = .init(dimMod.temperature, 1, .initUniqueSymbol("degF"));
pub const degR: Unit = defUnitFromUnit(si.K, 5 / 9, .initUniqueSymbol("degR"));
