const std = @import("std");
const testing = std.testing;
const si = @import("si.zig");
const Unit = @import("unit.zig").Unit;

// The SI base units, one per base dimension, indexed by dimension field name.
/// Used by `Unit.decompose(.SI)` to expand a unit into `m`, `kg`, `s`, etc.
pub const BASE_UNIT_SI: SystemUnit = .init(@constCast(&[8][]const u8{ "l", "m", "t", "i", "th", "n", "j", "a" }), @constCast(&[8]Unit{ si.m, si.kg, si.s, si.A, si.K, si.mol, si.cd, si.rad }));

// The CGS base units. Differs from SI in length (`cm`) and mass (`g`), and omits
/// electric current since CGS expresses electromagnetism through mechanical units.
pub const BASE_UNIT_CGS: SystemUnit = .init(@constCast(&[7][]const u8{ "l", "m", "t", "th", "n", "j", "a" }), @constCast(&[7]Unit{ si.cm, si.g, si.s, si.K, si.mol, si.cd, si.rad }));

/// Upper bound on the number of base dimensions a system can define (the eight
/// dimensions tracked by `Dim`). Backing arrays in `SystemUnit` are sized to this.
pub const MAX_SYSTEM_LENGTH = 8;

/// A fixed-capacity mapping from base-dimension field names (e.g. "l", "th") to the
/// base `Unit` representing that dimension in a particular system. Stored as two
/// parallel arrays of length `len` rather than a hash map, since the set is tiny and
/// known at comptime.
pub const SystemUnit = struct {
    const Self = @This();

    /// Number of populated entries in `dim_field` / `unit`.
    len: usize = MAX_SYSTEM_LENGTH,
    /// Dimension field names (parallel to `unit`).
    dim_field: [MAX_SYSTEM_LENGTH][]const u8 = undefined,
    /// Base units (parallel to `dim_field`).
    unit: [MAX_SYSTEM_LENGTH]Unit = undefined,

    // Builds a SystemUnit from parallel name/unit slices. Both slices must have the
    // same length and fit within MAX_SYSTEM_LENGTH; violations are a compile error at
    // comptime or a panic at runtime.
    fn init(dim_field: [][]const u8, unit: []Unit) Self {
        if (dim_field.len != unit.len) {
            if (@inComptime()) {
                @compileError("Trying to initialise BaseUnit with `dim_field` `and unit` of different lenght");
            } else {
                @panic("Trying to initialise BaseUnit with `dim_field` `and unit` of different lenght");
            }
        }
        if (dim_field.len > MAX_SYSTEM_LENGTH) {
            if (@inComptime()) {
                @compileError("`dim_field` length greater than 8 unsupported, panic !");
            } else {
                @panic("`dim_field` lenght greater than 8 unsupported, panic !");
            }
        }
        var result: Self = .{};
        result.len = dim_field.len;
        for (0..result.len) |i| {
            result.dim_field[i] = dim_field[i];
            result.unit[i] = unit[i];
        }
        return result;
    }

    // Returns the array index of a dimension field name, or null if absent. Field
    // names are at most two characters, so longer inputs short-circuit to null.
    fn getIdx(self: Self, field: []const u8) ?usize {
        if (field.len > 2) return null;
        for (0..self.len) |i| {
            if (std.mem.eql(u8, field, self.dim_field[i])) return i;
        }
        return null;
    }

    /// Returns the base unit for the given dimension field name, or null if the
    /// system does not define that dimension (e.g. "i" in CGS).
    pub fn get(self: Self, field: []const u8) ?Unit {
        if (self.getIdx(field)) |idx| {
            return self.unit[idx];
        }
        return null;
    }
};

/// The supported unit systems. `SI` and `CGS` have defined base units; the
/// electromagnetic variants are placeholders pending implementation.
pub const System = enum {
    const Self = @This();

    SI,
    CGS,
    ESU,
    EMU,
    GAUSS,

    /// Returns the base-unit table for this system. Panics for systems whose base
    /// units are not yet defined.
    pub fn get_base_unit(self: Self) SystemUnit {
        switch (self) {
            .SI => return BASE_UNIT_SI,
            .CGS => return BASE_UNIT_CGS,
            else => @panic("system has no base units."), // TODO: Make this a proper error !
        }
    }
};

test "SystemUnit init" {
    const si_unit = SystemUnit.init(
        @constCast(&[MAX_SYSTEM_LENGTH][]const u8{ "l", "m", "t", "i", "th", "n", "j", "a" }),
        @constCast(&[MAX_SYSTEM_LENGTH]Unit{ si.m, si.kg, si.s, si.A, si.K, si.mol, si.cd, si.rad }),
    );
    try testing.expectEqual(MAX_SYSTEM_LENGTH, si_unit.len);
    try testing.expectEqualStrings("i", si_unit.dim_field[3]);
    try testing.expectEqual(si.kg, si_unit.unit[1]);
}

test "SystemUnit getIdx" {
    const idx = BASE_UNIT_SI.getIdx("l").?;
    try testing.expectEqual(0, idx);
    const idx2 = BASE_UNIT_SI.getIdx("csa");
    try testing.expectEqual(null, idx2);
}

test "SystemUnit get" {
    const u_1: Unit = BASE_UNIT_SI.get("th").?;
    try testing.expectEqual(si.K, u_1);
    const u_2: ?Unit = BASE_UNIT_CGS.get("dfa");
    try testing.expectEqual(null, u_2);
}

test "System get_base_unit" {
    const sys: System = .SI;
    try testing.expectEqual(BASE_UNIT_SI, sys.get_base_unit());
    const sys_cgs: System = .CGS;
    try testing.expectEqual(BASE_UNIT_CGS, sys_cgs.get_base_unit());
}
