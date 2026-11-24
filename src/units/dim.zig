const std = @import("std");
const testing = std.testing;
const mem = std.mem;

const fraction = @import("fraction.zig");
const Fraction = fraction.Fraction;
const FractionError = fraction.FractionError;

// The eight base dimensions, each a `Dim` with a single exponent of 1. These are
// the building blocks from which all unit dimensions are composed. The first seven
// are the SI base dimensions; `angle` is included separately to support angular
// equivalencies even though angle is formally dimensionless in SI.
pub const length: Dim = .initUniqueFieldInt("l", 1);
pub const mass: Dim = .initUniqueFieldInt("m", 1);
pub const time: Dim = .initUniqueFieldInt("t", 1);
pub const electric_current: Dim = .initUniqueFieldInt("i", 1);
pub const temperature: Dim = .initUniqueFieldInt("th", 1);
pub const amount: Dim = .initUniqueFieldInt("n", 1);
pub const luminous_intensity: Dim = .initUniqueFieldInt("j", 1);
pub const angle: Dim = .initUniqueFieldInt("a", 1);

pub const ExponentInt = i32;
pub const Exponent = Fraction(ExponentInt);

/// The dimensional signature of a unit: the exponent of each base dimension,
/// stored as a `Fraction` so that fractional powers (e.g. from square roots) are
/// represented exactly. Two units are dimensionally compatible iff their `Dim`s
/// are equal.
///
/// Fields, each a `Exponent` exponent:
/// - `l`:  length
/// - `m`:  mass
/// - `t`:  time
/// - `i`:  electric current
/// - `th`: thermodynamic temperature
/// - `n`: amount of substance
/// - `j`:  luminous intensity
/// - `a`:  angle (dimensionless in SI, tracked here for angular equivalencies)
///
/// Arithmetic on dimensions mirrors arithmetic on units: multiplying two units adds
/// their dimensions, dividing subtracts them, raising to a power scales them. Each
/// operation returns `FractionError` if the underlying fraction arithmetic fails,
/// and comes in both a pure and an `InPlace` form.
pub const Dim = struct {
    const Self = @This();

    l: Exponent,
    m: Exponent,
    t: Exponent,
    i: Exponent,
    th: Exponent,
    n: Exponent,
    j: Exponent,
    a: Exponent,

    /// Constructs a dimensionless signature: every base-dimension exponent is zero.
    pub fn initDimensionless() Self {
        return .{
            .l = Exponent.initInt(0),
            .m = Exponent.initInt(0),
            .t = Exponent.initInt(0),
            .i = Exponent.initInt(0),
            .th = Exponent.initInt(0),
            .n = Exponent.initInt(0),
            .j = Exponent.initInt(0),
            .a = Exponent.initInt(0),
        };
    }

    /// Constructs a signature with a single base dimension (named by `field_name`,
    /// e.g. "l" or "th") set to `frac` and all others zero.
    pub fn initUniqueField(comptime field_name: []const u8, frac: Exponent) Self {
        var res = Self.initDimensionless();
        @field(res, field_name) = frac;
        return res;
    }

    /// Integer-exponent convenience wrapper around `initUniqueField`.
    pub fn initUniqueFieldInt(comptime field_name: []const u8, value: ExponentInt) Self {
        return Self.initUniqueField(field_name, Exponent.initInt(value));
    }

    /// Returns true if every base-dimension exponent matches, i.e. the two
    /// signatures are dimensionally identical.
    pub fn eql(self: Self, other: Self) bool {
        return (self.l.eql(other.l)) and
            (self.m.eql(other.m)) and
            (self.t.eql(other.t)) and
            (self.i.eql(other.i)) and
            (self.th.eql(other.th)) and
            (self.n.eql(other.n)) and
            (self.j.eql(other.j)) and
            (self.a.eql(other.a));
    }

    pub fn isDimensionless(self: Self) bool {
        return (self.l.isZero()) and
            (self.m.isZero()) and
            (self.t.isZero()) and
            (self.i.isZero()) and
            (self.th.isZero()) and
            (self.n.isZero()) and
            (self.j.isZero()) and
            (self.a.isZero());
    }

    /// Adds two signatures component-wise (the dimension of a unit product).
    /// Returns `FractionError` if any exponent addition overflows.
    pub fn add(self: Self, other: Self) FractionError!Self {
        var res: Self = undefined;
        inline for (@typeInfo(Self).@"struct".fields) |field| {
            const name = field.name;
            @field(res, name) = try @field(self, name).add(@field(other, name));
        }
        return res;
    }

    /// In-place variant of `add`.
    pub fn addInPlace(self: *Self, other: Self) FractionError!void {
        self.* = try self.add(other);
    }

    /// Subtracts two signatures component-wise (the dimension of a unit quotient).
    /// Returns `FractionError` if any exponent subtraction overflows.
    pub fn sub(self: Self, other: Self) FractionError!Self {
        var res: Self = undefined;
        inline for (@typeInfo(Self).@"struct".fields) |field| {
            const name = field.name;
            @field(res, name) = try @field(self, name).sub(@field(other, name));
        }
        return res;
    }

    /// In-place variant of `sub`.
    pub fn subInPlace(self: *Self, other: Self) FractionError!void {
        self.* = try self.sub(other);
    }

    /// Negates every exponent (the dimension of a reciprocal unit). Returns
    /// `FractionError` on fraction failure.
    pub fn neg(self: Self) FractionError!Self {
        var res = self;
        inline for (@typeInfo(Self).@"struct".fields) |field| {
            const name = field.name;
            @field(res, name) = try @field(self, name).neg();
        }
        return res;
    }

    /// In-place variant of `neg`.
    pub fn negInPlace(self: *Self) FractionError!void {
        self.* = try self.neg();
    }

    /// Multiplies every exponent by an integer (the dimension of a unit raised to
    /// an integer power). Returns `FractionError` on overflow.
    pub fn mulScalar(self: Self, value: ExponentInt) FractionError!Self {
        var res = self;
        inline for (@typeInfo(Self).@"struct".fields) |field| {
            const name = field.name;
            @field(res, name) = try @field(res, name).mulScalar(value);
        }
        return res;
    }

    /// In-place variant of `mulScalar`.
    pub fn mulScalarInPlace(self: *Self, value: ExponentInt) FractionError!void {
        self.* = try self.mulScalar(value);
    }

    /// Multiplies every exponent by a fraction (the dimension of a unit raised to a
    /// rational power). Returns `FractionError` on overflow.
    pub fn mulFraction(self: Self, frac: Exponent) FractionError!Self {
        var res = self;
        inline for (@typeInfo(Self).@"struct".fields) |field| {
            const name = field.name;
            @field(res, name) = try @field(res, name).mul(frac);
        }
        return res;
    }

    /// In-place variant of `mulFraction`.
    pub fn mulFractionInPlace(self: *Self, frac: Exponent) FractionError!void {
        self.* = try self.mulFraction(frac);
    }

    /// Divides every exponent by an integer (the dimension of an integer root, e.g.
    /// dividing by 2 for a square root). Returns `FractionError` on failure.
    pub fn divScalar(self: Self, value: ExponentInt) FractionError!Self {
        var res = self;
        inline for (@typeInfo(Self).@"struct".fields) |field| {
            const name = field.name;
            @field(res, name) = try @field(self, name).divScalar(value);
        }
        return res;
    }

    /// In-place variant of `divScalar`.
    pub fn divScalarInPlace(self: *Self, value: ExponentInt) FractionError!void {
        self.* = try self.divScalar(value);
    }

    /// Divides every exponent by a fraction. Returns `FractionError` on failure.
    pub fn divFraction(self: Self, frac: Exponent) FractionError!Self {
        var res = self;
        inline for (@typeInfo(Self).@"struct".fields) |field| {
            const name = field.name;
            @field(res, name) = try @field(res, name).div(frac);
        }
        return res;
    }

    /// In-place variant of `divFraction`.
    pub fn divFractionInPlace(self: *Self, frac: Exponent) FractionError!void {
        self.* = try self.divFraction(frac);
    }

    /// `std.fmt` formatting hook. Prints all eight exponents in a fixed order, e.g.
    /// `Dim: {l = 1, m = 0, t = 1, i = 0, th = 0, n = 0, j = 0, a = 0}`.
    pub fn format(self: Self, writer: *std.Io.Writer) std.Io.Writer.Error!void {
        try writer.writeAll("Dim: {");
        const fields = @typeInfo(Self).@"struct".fields;
        inline for (fields, 0..fields.len) |field, i| {
            const name = field.name;
            try writer.print("{s} = {f}", .{ name, @field(self, name) });
            if (i < fields.len - 1) {
                try writer.writeAll(", ");
            }
        }
        try writer.writeByte('}');
    }
};

test "init dimensionless" {
    const d = Dim.initDimensionless();
    inline for (@typeInfo(Dim).@"struct".fields) |field| {
        const name = field.name;
        try testing.expectEqual(Exponent.initInt(0), @field(d, name));
    }
}

test "initUniqueField" {
    const fieldValue = try Exponent.init(3, 4);
    const d = Dim.initUniqueField("t", fieldValue);
    inline for (@typeInfo(Dim).@"struct".fields) |field| {
        const name = field.name;
        if (mem.eql(u8, name, "t")) {
            try testing.expectEqual(try Exponent.init(3, 4), @field(d, name));
        } else {
            try testing.expectEqual(Exponent.initInt(0), @field(d, name));
        }
    }
}

test "initUniqueFieldInt" {
    const d = Dim.initUniqueFieldInt("l", 3);
    inline for (@typeInfo(Dim).@"struct".fields) |field| {
        const name = field.name;
        if (mem.eql(u8, name, "l")) {
            try testing.expectEqual(try Exponent.init(3, 1), @field(d, name));
        } else {
            try testing.expectEqual(Exponent.initInt(0), @field(d, name));
        }
    }
}

test "eql" {
    const d1 = Dim.initDimensionless();
    const d2 = Dim.initDimensionless();
    try testing.expect(Dim.eql(d1, d2) == true);
}

test "add" {
    var lengthTime: Dim = .initDimensionless();
    lengthTime.l = Exponent.initInt(1);
    lengthTime.t = Exponent.initInt(1);
    const added = try length.add(time);
    try testing.expectEqual(lengthTime, added);
}

test "addInPlace" {
    var massElectricCurrent: Dim = .initDimensionless();
    massElectricCurrent.m = Exponent.initInt(1);
    massElectricCurrent.i = Exponent.initInt(1);
    var massTest: Dim = mass;
    try massTest.addInPlace(electric_current);
    try testing.expectEqual(massElectricCurrent, massTest);
}

test "sub" {
    var tempByAmount: Dim = .initDimensionless();
    tempByAmount.th = Exponent.initInt(1);
    tempByAmount.n = Exponent.initInt(-1);
    const subbed = try temperature.sub(amount);
    try testing.expectEqual(tempByAmount, subbed);
}

test "subInPlace" {
    var lumBylength: Dim = .initDimensionless();
    lumBylength.l = Exponent.initInt(-1);
    lumBylength.j = Exponent.initInt(1);
    var lumTest: Dim = luminous_intensity;
    try lumTest.subInPlace(length);
    try testing.expectEqual(lumBylength, lumTest);
}

test "neg" {
    const invMass = try mass.neg();
    const invMassRes: Dim = .initUniqueFieldInt("m", -1);
    try testing.expectEqual(invMassRes, invMass);
}

test "mulScalar" {
    const massTestmul2 = try mass.mulScalar(2);
    try testing.expect(massTestmul2.m.eqlScalar(2) == true);

    const massTestmul_2 = try mass.mulScalar(-2);
    try testing.expect(massTestmul_2.m.eqlScalar(-2) == true);
    try testing.expect(massTestmul_2.l.eqlScalar(0) == true);
}

test "mulScalarInPlace" {
    var massTest = mass;
    try massTest.mulScalarInPlace(2);
    try testing.expect(massTest.m.eqlScalar(2) == true);

    var massTest2 = mass;
    try massTest2.mulScalarInPlace(-2);
    try testing.expect(massTest2.m.eqlScalar(-2) == true);
    try testing.expect(massTest2.l.eqlScalar(0) == true);
}

test "mulFraction" {
    const sizeTest4 = try length.mulScalar(4);
    const sizeTest4MulFrac = try sizeTest4.mulFraction(try Exponent.init(1, 4));
    try testing.expectEqual(length, sizeTest4MulFrac);
}

test "mulFractionInPlace" {
    var sizeTest4 = try length.mulScalar(4);
    try sizeTest4.mulFractionInPlace(try Exponent.init(1, 4));
    try testing.expectEqual(length, sizeTest4);
}

test "divScalar" {
    const massTestmul2 = try mass.divScalar(2);
    try testing.expectEqual(try Exponent.init(1, 2), massTestmul2.m);

    const massTestmul_2 = try mass.divScalar(-2);
    try testing.expectEqual(try Exponent.init(-1, 2), massTestmul_2.m);
    try testing.expect(massTestmul_2.l.eqlScalar(0) == true);
}

test "divScalarInPlace" {
    var massTest = mass;
    try massTest.divScalarInPlace(2);
    try testing.expectEqual(try Exponent.init(1, 2), massTest.m);

    var massTest2 = mass;
    try massTest2.divScalarInPlace(-2);
    try testing.expectEqual(try Exponent.init(-1, 2), massTest2.m);
    try testing.expect(massTest2.l.eqlScalar(0) == true);
}

test "divFraction" {
    const sizeTest4 = try length.divScalar(4);
    const sizeTest4MulFrac = try sizeTest4.divFraction(try Exponent.init(1, 4));
    try testing.expectEqual(length, sizeTest4MulFrac);
}

test "divFractionInPlace" {
    var sizeTest4 = try length.divScalar(4);
    try sizeTest4.divFractionInPlace(try Exponent.init(1, 4));
    try testing.expectEqual(length, sizeTest4);
}

test "negInPlace" {
    var amountTest = amount;
    try amountTest.negInPlace();
    const invAmountRes: Dim = .initUniqueFieldInt("n", -1);
    try testing.expectEqual(invAmountRes, amountTest);
}

test "format" {
    var lengthTime: Dim = .initDimensionless();
    lengthTime.l = Exponent.initInt(1);
    lengthTime.t = Exponent.initInt(1);

    var buf: [64]u8 = undefined;
    const print_res = try std.fmt.bufPrint(&buf, "{f}", .{lengthTime});
    const expected: []const u8 = "Dim: {l = 1, m = 0, t = 1, i = 0, th = 0, n = 0, j = 0, a = 0}";
    try testing.expectEqualSlices(u8, expected, print_res);
}
