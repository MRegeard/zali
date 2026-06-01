const std = @import("std");
const testing = std.testing;
const math = std.math;
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
pub const electricCurrent: Dim = .initUniqueFieldInt("i", 1);
pub const temperature: Dim = .initUniqueFieldInt("th", 1);
pub const amount: Dim = .initUniqueFieldInt("n", 1);
pub const luminousIntensity: Dim = .initUniqueFieldInt("j", 1);
pub const angle: Dim = .initUniqueFieldInt("a", 1);

/// The dimensional signature of a unit: the exponent of each base dimension,
/// stored as a `Fraction` so that fractional powers (e.g. from square roots) are
/// represented exactly. Two units are dimensionally compatible iff their `Dim`s
/// are equal.
///
/// Fields, each a `Fraction(isize)` exponent:
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

    l: Fraction(isize),
    m: Fraction(isize),
    t: Fraction(isize),
    i: Fraction(isize),
    th: Fraction(isize),
    n: Fraction(isize),
    j: Fraction(isize),
    a: Fraction(isize),

    /// Constructs a dimensionless signature: every base-dimension exponent is zero.
    pub fn initDimensionless() Self {
        return .{
            .l = Fraction(isize).initInt(0),
            .m = Fraction(isize).initInt(0),
            .t = Fraction(isize).initInt(0),
            .i = Fraction(isize).initInt(0),
            .th = Fraction(isize).initInt(0),
            .n = Fraction(isize).initInt(0),
            .j = Fraction(isize).initInt(0),
            .a = Fraction(isize).initInt(0),
        };
    }

    /// Constructs a signature with a single base dimension (named by `field_name`,
    /// e.g. "l" or "th") set to `frac` and all others zero.
    pub fn initUniqueField(comptime field_name: []const u8, frac: Fraction(isize)) Self {
        var res = Self.initDimensionless();
        @field(res, field_name) = frac;
        return res;
    }

    /// Integer-exponent convenience wrapper around `initUniqueField`.
    pub fn initUniqueFieldInt(comptime field_name: []const u8, value: isize) Self {
        return Self.initUniqueField(field_name, Fraction(isize).initInt(value));
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

    /// Adds two signatures component-wise (the dimension of a unit product).
    /// Returns `FractionError` if any exponent addition overflows.
    pub fn add(self: Self, other: Self) FractionError!Self {
        var res = Self.initDimensionless();
        inline for (@typeInfo(Self).@"struct".fields) |field| {
            const name = field.name;
            @field(res, name) = try @field(self, name).add(@field(other, name));
        }
        return res;
    }

    /// In-place variant of `add`.
    pub fn addInPlace(self: *Self, other: Self) FractionError!void {
        inline for (@typeInfo(Self).@"struct".fields) |field| {
            const name = field.name;
            try @field(self.*, name).addInPlace(@field(other, name));
        }
    }

    /// Subtracts two signatures component-wise (the dimension of a unit quotient).
    /// Returns `FractionError` if any exponent subtraction overflows.
    pub fn sub(self: Self, other: Self) FractionError!Self {
        var res = Self.initDimensionless();
        inline for (@typeInfo(Self).@"struct".fields) |field| {
            const name = field.name;
            @field(res, name) = try @field(self, name).sub(@field(other, name));
        }
        return res;
    }

    /// In-place variant of `sub`.
    pub fn subInPlace(self: *Self, other: Self) FractionError!void {
        inline for (@typeInfo(Self).@"struct".fields) |field| {
            const name = field.name;
            try @field(self.*, name).subInPlace(@field(other, name));
        }
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
        inline for (@typeInfo(Self).@"struct".fields) |field| {
            const name = field.name;
            try @field(self.*, name).negInPlace();
        }
    }

    /// Multiplies every exponent by an integer (the dimension of a unit raised to
    /// an integer power). Returns `FractionError` on overflow.
    pub fn mulScalar(self: Self, value: isize) FractionError!Self {
        var res = self;
        inline for (@typeInfo(Self).@"struct".fields) |field| {
            const name = field.name;
            @field(res, name) = try @field(res, name).mulScalar(value);
        }
        return res;
    }

    /// In-place variant of `mulScalar`.
    pub fn mulScalarInPlace(self: *Self, value: isize) FractionError!void {
        inline for (@typeInfo(Self).@"struct".fields) |field| {
            const name = field.name;
            try @field(self.*, name).mulScalarInPlace(value);
        }
    }

    /// Multiplies every exponent by a fraction (the dimension of a unit raised to a
    /// rational power). Returns `FractionError` on overflow.
    pub fn mulFraction(self: Self, frac: Fraction(isize)) FractionError!Self {
        var res = self;
        inline for (@typeInfo(Self).@"struct".fields) |field| {
            const name = field.name;
            @field(res, name) = try @field(res, name).mul(frac);
        }
        return res;
    }

    /// In-place variant of `mulFraction`.
    pub fn mulFractionInPlace(self: *Self, frac: Fraction(isize)) FractionError!void {
        inline for (@typeInfo(Self).@"struct".fields) |field| {
            const name = field.name;
            try @field(self.*, name).mulInPlace(frac);
        }
    }

    /// Divides every exponent by an integer (the dimension of an integer root, e.g.
    /// dividing by 2 for a square root). Returns `FractionError` on failure.
    pub fn divScalar(self: Self, value: isize) FractionError!Self {
        var res = self;
        inline for (@typeInfo(Self).@"struct".fields) |field| {
            const name = field.name;
            @field(res, name) = try @field(self, name).divScalar(value);
        }
        return res;
    }

    /// In-place variant of `divScalar`.
    pub fn divScalarInPlace(self: *Self, value: isize) FractionError!void {
        inline for (@typeInfo(Self).@"struct".fields) |field| {
            const name = field.name;
            try @field(self.*, name).divScalarInPlace(value);
        }
    }

    /// Divides every exponent by a fraction. Returns `FractionError` on failure.
    pub fn divFraction(self: Self, frac: Fraction(isize)) FractionError!Self {
        var res = self;
        inline for (@typeInfo(Self).@"struct".fields) |field| {
            const name = field.name;
            @field(res, name) = try @field(res, name).div(frac);
        }
        return res;
    }

    /// In-place variant of `divFraction`.
    pub fn divFractionInPlace(self: *Self, frac: Fraction(isize)) FractionError!void {
        inline for (@typeInfo(Self).@"struct".fields) |field| {
            const name = field.name;
            try @field(self.*, name).divInPlace(frac);
        }
    }

    /// `std.fmt` formatting hook. Prints all eight exponents in a fixed order, e.g.
    /// `Dim: {l = 1, m = 0, t = 1, i = 0, th = 0, n = 0, j = 0, a = 0}`.
    pub fn format(self: *const Self, writer: *std.Io.Writer) std.Io.Writer.Error!void {
        try writer.writeAll("Dim: {");
        inline for (@typeInfo(Self).@"struct".fields, 0..8) |field, i| {
            const name = field.name;
            try writer.print("{s} = {f}", .{ name, @field(self.*, name) });
            if (i < 7) {
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
        try testing.expectEqual(Fraction(isize).initInt(0), @field(d, name));
    }
}

test "initUniqueField" {
    const fieldValue = try Fraction(isize).init(3, 4);
    const d = Dim.initUniqueField("t", fieldValue);
    inline for (@typeInfo(Dim).@"struct".fields) |field| {
        const name = field.name;
        if (mem.eql(u8, name, "t")) {
            try testing.expectEqual(try Fraction(isize).init(3, 4), @field(d, name));
        } else {
            try testing.expectEqual(Fraction(isize).initInt(0), @field(d, name));
        }
    }
}

test "initUniqueFieldInt" {
    const d = Dim.initUniqueFieldInt("l", 3);
    inline for (@typeInfo(Dim).@"struct".fields) |field| {
        const name = field.name;
        if (mem.eql(u8, name, "l")) {
            try testing.expectEqual(try Fraction(isize).init(3, 1), @field(d, name));
        } else {
            try testing.expectEqual(Fraction(isize).initInt(0), @field(d, name));
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
    lengthTime.l = Fraction(isize).initInt(1);
    lengthTime.t = Fraction(isize).initInt(1);
    const added = try length.add(time);
    try testing.expectEqual(lengthTime, added);
}

test "addInPlace" {
    var massElectricCurrent: Dim = .initDimensionless();
    massElectricCurrent.m = Fraction(isize).initInt(1);
    massElectricCurrent.i = Fraction(isize).initInt(1);
    var massTest: Dim = mass;
    try massTest.addInPlace(electricCurrent);
    try testing.expectEqual(massElectricCurrent, massTest);
}

test "sub" {
    var tempByAmount: Dim = .initDimensionless();
    tempByAmount.th = Fraction(isize).initInt(1);
    tempByAmount.n = Fraction(isize).initInt(-1);
    const subbed = try temperature.sub(amount);
    try testing.expectEqual(tempByAmount, subbed);
}

test "subInPlace" {
    var lumBylength: Dim = .initDimensionless();
    lumBylength.l = Fraction(isize).initInt(-1);
    lumBylength.j = Fraction(isize).initInt(1);
    var lumTest: Dim = luminousIntensity;
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
    const sizeTest4MulFrac = try sizeTest4.mulFraction(try Fraction(isize).init(1, 4));
    try testing.expectEqual(length, sizeTest4MulFrac);
}

test "mulFractionInPlace" {
    var sizeTest4 = try length.mulScalar(4);
    try sizeTest4.mulFractionInPlace(try Fraction(isize).init(1, 4));
    try testing.expectEqual(length, sizeTest4);
}

test "divScalar" {
    const massTestmul2 = try mass.divScalar(2);
    try testing.expectEqual(try Fraction(isize).init(1, 2), massTestmul2.m);

    const massTestmul_2 = try mass.divScalar(-2);
    try testing.expectEqual(try Fraction(isize).init(-1, 2), massTestmul_2.m);
    try testing.expect(massTestmul_2.l.eqlScalar(0) == true);
}

test "divScalarInPlace" {
    var massTest = mass;
    try massTest.divScalarInPlace(2);
    try testing.expectEqual(try Fraction(isize).init(1, 2), massTest.m);

    var massTest2 = mass;
    try massTest2.divScalarInPlace(-2);
    try testing.expectEqual(try Fraction(isize).init(-1, 2), massTest2.m);
    try testing.expect(massTest2.l.eqlScalar(0) == true);
}

test "divFraction" {
    const sizeTest4 = try length.divScalar(4);
    const sizeTest4MulFrac = try sizeTest4.divFraction(try Fraction(isize).init(1, 4));
    try testing.expectEqual(length, sizeTest4MulFrac);
}

test "divFractionInPlace" {
    var sizeTest4 = try length.divScalar(4);
    try sizeTest4.divFractionInPlace(try Fraction(isize).init(1, 4));
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
    lengthTime.l = Fraction(isize).initInt(1);
    lengthTime.t = Fraction(isize).initInt(1);

    var buf: [64]u8 = undefined;
    const print_res = try std.fmt.bufPrint(&buf, "{f}", .{lengthTime});
    const expected: []const u8 = "Dim: {l = 1, m = 0, t = 1, i = 0, th = 0, n = 0, j = 0, a = 0}";
    try testing.expectEqualSlices(u8, expected, print_res);
}
