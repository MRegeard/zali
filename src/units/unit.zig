const std = @import("std");
const fmt = std.fmt;
const math = std.math;
const mem = std.mem;
const testing = std.testing;
const dimMod = @import("dim.zig");
const Dim = dimMod.Dim;
const fraction = @import("fraction.zig");
const symbolMod = @import("symbol.zig");
const SymbolExpression = symbolMod.SymbolExpression;
const SymbolTerm = symbolMod.SymbolTerm;
const systemMod = @import("system.zig");
const System = systemMod.System;
const SystemUnit = systemMod.SystemUnit;
const utils = @import("utils.zig");

const FractionError = fraction.FractionError;
const Fraction = fraction.Fraction;

/// The dimensionless unit: zero in every physical dimension, scale 1, empty symbol.
/// Used as the identity for unit composition and as the unit of pure scalars.
pub const UNITLESS: Unit = .init(Dim.initDimensionless(), 1.0, .initUniqueSymbol(""));

/// A unit of measurement. Carries everything needed to
/// describe and convert a physical unit:
/// - `dim`:    the dimensional signature (length, mass, time, ... each as a `Fraction`
///             exponent), used for compatibility checks and dimensional analysis.
/// - `scale`:  the multiplicative factor relative to the base unit of the dimension
///             (e.g. `m` -> 1, `cm` -> 0.01).
/// - `offset`: an optional affine offset for non-multiplicative scales such as
///             Celsius or Fahrenheit (`null` for ordinary multiplicative units).
/// - `symbol`: the symbolic representation, a product of named terms with exponents.
///
/// Most composition operations come in two forms: a `try*` variant that returns a
/// `FractionError` when the underlying dimension arithmetic could overflow, and a
/// plain variant that treats such failure as `unreachable` (safe because physical
/// exponents are small). Affine units cannot be composed; attempting to do so is a
/// compile error at comptime or a panic at runtime.
pub const Unit = struct {
    const Self = @This();

    dim: Dim,
    scale: f64 = 1.0,
    offset: ?f64 = null,
    symbol: SymbolExpression,

    /// Constructs an ordinary (non-affine) unit with no offset.
    pub fn init(dim: Dim, scale: f64, symbol: SymbolExpression) Self {
        return .{ .dim = dim, .scale = scale, .symbol = symbol, .offset = null };
    }

    /// Constructs an affine unit with an offset (e.g. Celsius, Fahrenheit), where
    /// conversion is `value * scale + offset` rather than a pure multiplication.
    pub fn initAffine(dim: Dim, scale: f64, offset: f64, symbol: SymbolExpression) Self {
        return .{ .dim = dim, .scale = scale, .offset = offset, .symbol = symbol };
    }

    /// Strict equality: compares dimension, scale, offset, AND symbol. Two units that
    /// are physically identical but spelled differently are NOT equal under this.
    pub fn eqlExact(self: Self, other: Self) bool {
        return self.dim.eql(other.dim) and
            self.scale == other.scale and
            self.offset == other.offset and
            self.symbol.eql(other.symbol);
    }

    /// Physical equality: compares dimension, scale, and offset but ignores the
    /// symbol. Units that represent the same physical quantity compare equal even if
    /// their symbolic spelling differs.
    pub fn eql(self: Self, other: Self) bool {
        return self.dim.eql(other.dim) and self.scale == other.scale and self.offset == other.offset;
    }

    /// Multiplies two units: adds their dimensions, multiplies their scales, and
    /// concatenates their symbols. Returns `FractionError` if the dimension addition
    /// overflows. Affine units (non-null offset) cannot be multiplied.
    pub fn tryMul(self: Self, other: Self) FractionError!Self {
        if (self.offset != null or other.offset != null) {
            if (@inComptime()) {
                @compileError("Cannot multiply affine units (non-zero offset).");
            } else {
                @panic("Cannot multiply affine units (non-zero offset).");
            }
        }
        return .{
            .dim = try self.dim.add(other.dim),
            .scale = self.scale * other.scale,
            .offset = null,
            .symbol = self.symbol.concatenate(other.symbol),
        };
    }

    /// Infallible `tryMul`: treats dimension overflow as `unreachable`.
    pub fn mul(self: Self, other: Self) Self {
        const res = self.tryMul(other) catch unreachable;
        return res;
    }

    /// Divides two units: subtracts their dimensions, divides their scales, and
    /// concatenates `self`'s symbol with the inverse of `other`'s. Returns
    /// `FractionError` on dimension overflow. Affine units cannot be divided.
    pub fn tryDiv(self: Self, other: Self) FractionError!Self {
        if (self.offset != null or other.offset != null) {
            if (@inComptime()) {
                @compileError("Cannot divide affine units (non-zero offset).");
            } else {
                @panic("Cannot divide affine units (non-zero offset).");
            }
        }
        return .{
            .dim = try self.dim.sub(other.dim),
            .scale = self.scale / other.scale,
            .offset = null,
            .symbol = self.symbol.concatenate(other.symbol.inv()),
        };
    }

    /// Infallible `tryDiv`: treats dimension overflow as `unreachable`.
    pub fn div(self: Self, other: Self) Self {
        return self.tryDiv(other) catch unreachable;
    }

    /// Raises a unit to an integer power: multiplies each dimension exponent by
    /// `value`, raises the scale to that power, and scales the symbol exponents.
    /// Returns `FractionError` on dimension overflow. Affine units cannot be powered.
    pub fn tryPow(self: Self, value: isize) FractionError!Self {
        if (self.offset != null) {
            if (@inComptime()) {
                @compileError("Cannot power affine units (non-zero offset).");
            } else {
                @panic("Cannot power affine units (non-zero offset).");
            }
        }
        const value_f: f64 = @floatFromInt(value);
        return .{
            .dim = try self.dim.mulScalar(value),
            .scale = math.pow(f64, self.scale, value_f),
            .offset = null,
            .symbol = self.symbol.mulFrac(.initInt(value)),
        };
    }

    /// Infallible `tryPow`: treats dimension overflow as `unreachable`.
    pub fn pow(self: Self, value: isize) Self {
        return self.tryPow(value) catch unreachable;
    }

    /// Square root of a unit: halves each dimension exponent, takes the square root
    /// of the scale, and halves the symbol exponents. Returns `FractionError` if a
    /// dimension exponent isn't evenly halvable in the fraction representation.
    /// Affine units cannot be rooted.
    pub fn trySqrt(self: Self) FractionError!Self {
        if (self.offset != null) {
            if (@inComptime()) {
                @compileError("Cannot square root affine units (non-zero offset).");
            } else {
                @panic("Cannot square root affine units (non-zero offset).");
            }
        }
        return .{
            .dim = try self.dim.divScalar(2),
            .scale = math.sqrt(self.scale),
            .offset = null,
            .symbol = self.symbol.divFrac(.initInt(2)),
        };
    }

    /// Infallible `trySqrt`: treats fraction failure as `unreachable`.
    pub fn sqrt(self: Self) Self {
        return self.trySqrt() catch unreachable;
    }

    /// Cube root of a unit: divides each dimension exponent by three, takes the cube
    /// root of the scale, and divides the symbol exponents by three. Returns
    /// `FractionError` on fraction failure. Affine units cannot be rooted.
    pub fn tryCbrt(self: Self) FractionError!Self {
        if (self.offset != null) {
            if (@inComptime()) {
                @compileError("Cannot cube root affine units (non-zero offset).");
            } else {
                @panic("Cannot cube root affine units (non-zero offset).");
            }
        }
        return .{
            .dim = try self.dim.divScalar(3),
            .scale = math.cbrt(self.scale),
            .offset = null,
            .symbol = self.symbol.divFrac(.initInt(3)),
        };
    }

    /// Infallible `tryCbrt`: treats fraction failure as `unreachable`.
    pub fn cbrt(self: Self) Self {
        return self.tryCbrt() catch unreachable;
    }

    /// Raises a unit to a rational power `num/denum`: multiplies each dimension
    /// exponent by the fraction, raises the scale to the fractional power, and scales
    /// the symbol exponents. Returns `FractionError` on invalid fraction or overflow.
    /// Affine units cannot be powered.
    pub fn tryPowByFraction(self: Self, num: isize, denum: isize) FractionError!Self {
        if (self.offset != null) {
            if (@inComptime()) {
                @compileError("Cannot power affine units (non-zero offset).");
            } else {
                @panic("Cannot power affine units (non-zero offset).");
            }
        }
        const frac = try Fraction(isize).init(num, denum);
        return .{
            .dim = try self.dim.mulFraction(frac),
            .scale = math.pow(f64, self.scale, frac.toFloat(f64)),
            .offset = null,
            .symbol = self.symbol.mulFrac(frac),
        };
    }

    /// Infallible `tryPowByFraction`: treats fraction failure as `unreachable`.
    pub fn powByFraction(self: Self, num: isize, denum: isize) Self {
        return self.tryPowByFraction(num, denum) catch unreachable;
    }

    /// Like `tryPowByFraction` but takes an already-constructed `Fraction` directly,
    /// avoiding re-validating the numerator/denominator. Returns `FractionError` on
    /// dimension overflow. Affine units cannot be powered.
    pub fn tryPowByZaliFraction(self: Self, frac: Fraction(isize)) FractionError!Self {
        if (self.offset != null) {
            if (@inComptime()) {
                @compileError("Cannot power affine units (non-zero offset).");
            } else {
                @panic("Cannot power affine units (non-zero offset).");
            }
        }
        return .{
            .dim = try self.dim.mulFraction(frac),
            .scale = math.pow(f64, self.scale, frac.toFloat(f64)),
            .offset = null,
            .symbol = self.symbol.mulFrac(frac),
        };
    }

    /// Infallible `tryPowByZaliFraction`: treats fraction failure as `unreachable`.
    pub fn powByZaliFraction(self: Self, frac: Fraction(isize)) Self {
        return self.tryPowByZaliFraction(frac) catch unreachable;
    }

    /// `std.fmt` formatting hook. Prints a multi-line description of the unit: each
    /// non-zero dimension exponent (as integer or fraction), then the scale, offset,
    /// and symbol.
    pub fn format(self: Self, writer: *std.Io.Writer) std.Io.Writer.Error!void {
        try writer.writeAll("Unit:\n    dim:");
        inline for (@typeInfo(Dim).@"struct".fields) |field| {
            const name = field.name;
            var field_value = @field(self.dim, name);
            field_value.reduceInPlace();
            if (!field_value.eqlScalar(0)) {
                if (field_value.isInt()) {
                    try writer.print(" {s} = {d},", .{ name, field_value.num });
                } else {
                    try writer.print(" {s} = {d}/{d},", .{ name, field_value.num, field_value.denum });
                }
            }
        }
        try writer.print("\n    scale: {d}\n    offset: {?d}\n    symbols: {f}", .{ self.scale, self.offset, self.symbol });
    }

    /// Decomposes the unit into the base units of `system` (e.g. N -> kg m s-2).
    /// The returned unit has a scale of 1. The original scale is applied to the
    /// value instead when decomposing a quantity, so 20 cm becomes 0.2 m.
    pub fn decompose(self: Self, system: System) Self {
        var new_unit: Self = UNITLESS;
        const system_unit = system.get_base_unit();
        inline for (@typeInfo(Dim).@"struct".fields) |field| {
            const name = field.name;
            const assoc_unit = system_unit.get(name) orelse {
                utils.compErrOrPanic("Cannot decompose unit into the requested system");
                unreachable;
            };
            const dim_field_value = @field(self.dim, name);
            if (!dim_field_value.eqlScalar(0)) {
                @field(new_unit.dim, name) = dim_field_value;
                if (assoc_unit.symbol.len != 1) utils.compErrOrPanic("Is this even possible ?");
                const new_unit_symbol: SymbolTerm = .init(assoc_unit.symbol.terms[0].symbol, dim_field_value);
                new_unit.symbol.appendInPlace(new_unit_symbol);
            }
        }
        return new_unit;
    }
};

test "test init" {
    const meter: Unit = .init(dimMod.length, 1.0, try SymbolExpression.initFromString("m"));
    try testing.expectEqual(dimMod.length, meter.dim);
    try testing.expectEqual(1.0, meter.scale);
    try testing.expectEqualSlices(u8, "m", meter.symbol.terms[0].symbol);
    try testing.expect(meter.symbol.terms[0].exponent.eqlScalar(1));
    try testing.expectEqual(null, meter.offset);
}

test "test init affine" {
    const degC: Unit = .initAffine(dimMod.temperature, 1.0, 273.15, try SymbolExpression.initFromString("degC"));
    try testing.expectEqual(dimMod.temperature, degC.dim);
    try testing.expectEqual(1.0, degC.scale);
    try testing.expectEqualSlices(u8, "degC", degC.symbol.terms[0].symbol);
    try testing.expectEqual(273.15, degC.offset);
}

test "test multiply units" {
    const meter1: Unit = .init(dimMod.length, 1, try SymbolExpression.initFromString("m"));
    const meter2: Unit = .init(dimMod.length, 1, try SymbolExpression.initFromString("m"));
    const meter_square = try meter1.tryMul(meter2);
    try testing.expectEqual(1.0, meter_square.scale);
    try testing.expectEqual(null, meter_square.offset);
    try testing.expectEqual(1, meter_square.symbol.len);
    try testing.expectEqualSlices(u8, "m", meter_square.symbol.terms[0].symbol);
    try testing.expect(meter_square.symbol.terms[0].exponent.eqlScalar(2));
    try testing.expectEqual(try dimMod.length.add(dimMod.length), meter_square.dim);
}

test "test divid units" {
    const meter: Unit = .init(dimMod.length, 1, try SymbolExpression.initFromString("m"));
    const centimeter: Unit = .init(dimMod.length, 0.01, try SymbolExpression.initFromString("cm"));
    const meter_centimeter = try meter.tryDiv(centimeter);
    try testing.expectEqual(100, meter_centimeter.scale);
    try testing.expectEqual(null, meter_centimeter.offset);
    try testing.expectEqual(2, meter_centimeter.symbol.len);
    try testing.expect(meter_centimeter.symbol.terms[1].exponent.eqlScalar(-1));
    try testing.expectEqualSlices(u8, "m", meter_centimeter.symbol.terms[0].symbol);
    try testing.expectEqual(try dimMod.length.sub(dimMod.length), meter_centimeter.dim);
}

test "pow" {
    const cm: Unit = .init(dimMod.length, 0.01, try SymbolExpression.initFromString("cm"));
    const cm_cube = try cm.tryPow(3);
    try testing.expect(math.approxEqAbs(f64, 1e-6, cm_cube.scale, 1e-15));
    try testing.expectEqual(null, cm_cube.offset);
    try testing.expectEqual(try dimMod.length.add(try dimMod.length.add(dimMod.length)), cm_cube.dim);
    try testing.expectEqual(1, cm_cube.symbol.len);
    try testing.expect(cm_cube.symbol.terms[0].exponent.eqlScalar(3));
    try testing.expectEqualSlices(u8, "cm", cm_cube.symbol.terms[0].symbol);
}

test "sqrt" {
    const m2: Unit = .init(try dimMod.length.add(dimMod.length), 1, try SymbolExpression.initFromString("m2"));
    const m: Unit = .init(dimMod.length, 1, try SymbolExpression.initFromString("m"));
    const m2sqrt = try m2.trySqrt();
    try testing.expect(math.approxEqAbs(f64, m.scale, m2sqrt.scale, 1e-15));
    try testing.expectEqual(null, m2sqrt.offset);
    try testing.expectEqual(m.dim, m2sqrt.dim);
    try testing.expectEqual(1, m2sqrt.symbol.len);
    try testing.expect(m2sqrt.symbol.terms[0].exponent.eqlScalar(1));
    try testing.expectEqualSlices(u8, "m", m2sqrt.symbol.terms[0].symbol);
}

test "cbrt" {
    const m3: Unit = .init(try dimMod.length.mulScalar(3), 1, try SymbolExpression.initFromString("m3"));
    const m: Unit = .init(dimMod.length, 1, try SymbolExpression.initFromString("m"));
    const m3cbrt = try m3.tryCbrt();
    try testing.expect(math.approxEqAbs(f64, m.scale, m3cbrt.scale, 1e-15));
    try testing.expectEqual(null, m3cbrt.offset);
    try testing.expectEqual(m.dim, m3cbrt.dim);
    try testing.expectEqual(1, m3cbrt.symbol.len);
    try testing.expect(m3cbrt.symbol.terms[0].exponent.eqlScalar(1));
    try testing.expectEqualSlices(u8, "m", m3cbrt.symbol.terms[0].symbol);
}

test "powByFraction" {
    const m4: Unit = .init(try dimMod.length.mulScalar(4), 1, try SymbolExpression.initFromString("m4"));
    const m: Unit = .init(dimMod.length, 1, try SymbolExpression.initFromString("m"));
    const m4PowByFrac = try m4.tryPowByFraction(1, 4);
    try testing.expect(math.approxEqAbs(f64, m.scale, m4PowByFrac.scale, 1e-15));
    try testing.expectEqual(null, m4PowByFrac.offset);
    try testing.expectEqual(m.dim, m4PowByFrac.dim);
    try testing.expectEqual(1, m4PowByFrac.symbol.len);
    try testing.expect(m4PowByFrac.symbol.terms[0].exponent.eqlScalar(1));
    try testing.expectEqualSlices(u8, "m", m4PowByFrac.symbol.terms[0].symbol);
}

test "powByZaliFraction" {
    const frac = try Fraction(isize).init(1, 4);
    const m4: Unit = .init(try dimMod.length.mulScalar(4), 1, try SymbolExpression.initFromString("m4"));
    const m: Unit = .init(dimMod.length, 1, try SymbolExpression.initFromString("m"));
    const m4PowByFrac = try m4.tryPowByZaliFraction(frac);
    try testing.expect(math.approxEqAbs(f64, m.scale, m4PowByFrac.scale, 1e-15));
    try testing.expectEqual(null, m4PowByFrac.offset);
    try testing.expectEqual(m.dim, m4PowByFrac.dim);
    try testing.expectEqual(1, m4PowByFrac.symbol.len);
    try testing.expect(m4PowByFrac.symbol.terms[0].exponent.eqlScalar(1));
    try testing.expectEqualSlices(u8, "m", m4PowByFrac.symbol.terms[0].symbol);
}

test "format" {
    const allocator = std.testing.allocator;
    const m4: Unit = .init(try dimMod.length.mulScalar(4), 1, try SymbolExpression.initFromString("m4"));
    const expected: []const u8 =
        \\Unit:
        \\    dim: l = 4,
        \\    scale: 1
        \\    offset: null
        \\    symbols: m4
    ;
    const print = try std.fmt.allocPrint(allocator, "{f}", .{m4});
    defer allocator.free(print);
    try testing.expectEqualStrings(expected, print);
}

test "decompose" {
    const dimension: Dim = .{
        .i = Fraction(isize).initInt(0),
        .j = Fraction(isize).initInt(0),
        .l = Fraction(isize).initInt(1),
        .m = Fraction(isize).initInt(1),
        .n = Fraction(isize).initInt(0),
        .t = Fraction(isize).initInt(-2),
        .th = Fraction(isize).initInt(0),
        .a = Fraction(isize).initInt(0),
    };
    const N: Unit = .init(dimension, 1, try SymbolExpression.initFromString("N"));
    const N_decomp = N.decompose(.SI);
    try testing.expectEqual(dimension, N_decomp.dim);
    try testing.expect((try SymbolExpression.initFromString("m kg s-2")).eql(N_decomp.symbol));
    try testing.expectApproxEqAbs(1.0, N_decomp.scale, 1e-15);
}
