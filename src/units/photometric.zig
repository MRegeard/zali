const std = @import("std");
const cst = @import("../constants.zig");
const unit = @import("unit.zig");
const Unit = unit.Unit;
const astrophys = @import("astrophys.zig");
const si = @import("si.zig");
const cgs = @import("cgs.zig");
const utils = @import("utils.zig");
const defUnitFromUnit = utils.defUnitFromUnit;
const defUnitFromConst = utils.defUnitFromConst;

const PI_F64: f64 = std.math.pi;

pub const Bol: Unit = defUnitFromConst(cst.L_bol0, .initUniqueSymbol("Bol"));
pub const bol: Unit = defUnitFromUnit(
    @TypeOf(cst.L_bol0.quantity).unit.div(astrophys.pc.pow(2)),
    cst.L_bol0.quantity.value / (4 * PI_F64 * 100),
    .initUniqueSymbol("bol"),
);
pub const AB: Unit = defUnitFromUnit(cgs.erg.mul(si.cm.pow(-2)).div(si.s).div(si.Hz), std.math.pow(f64, 10, 48.6 / -2.5), .initUniqueSymbol("AB"));
pub const ST: Unit = defUnitFromUnit(cgs.erg.mul(si.cm.pow(-2)).div(si.s).div(si.AA), std.math.pow(f64, 10, 21.1 / -2.5), .initUniqueSymbol("ST"));
pub const mgy: Unit = defUnitFromUnit(unit.UNITLESS, 1.0, .initUniqueSymbol("mgy"));
