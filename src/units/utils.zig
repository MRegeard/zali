const std = @import("std");
const constants = @import("../constants.zig");
const Constant = constants.Constant;
const Unit = @import("unit.zig").Unit;
const SymbolExpression = @import("symbol.zig").SymbolExpression;

pub fn defUnitFromUnit(
    comptime root: Unit,
    conversion_scale: f64,
    symbol: SymbolExpression,
) Unit {
    if (root.offset != null) {
        @compileError(
            "defUnitFromUnit() cannot be used with affine units (non-zero offset).",
        );
    }
    return Unit.init(root.dim, conversion_scale * root.scale, symbol);
}

pub fn defUnitFromConst(
    comptime constant: anytype,
    symbol: SymbolExpression,
) Unit {
    const unit: Unit = @TypeOf(constant.quantity).unit;
    return defUnitFromUnit(unit, constant.quantity.value, symbol);
}
