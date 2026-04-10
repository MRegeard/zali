const std = @import("std");
const testing = std.testing;
const si = @import("si.zig");
const Unit = @import("unit.zig").Unit;

pub const BASE_UNIT_SI: SystemUnit = .init(@constCast(&[8][]const u8{ "l", "m", "t", "i", "th", "n", "j", "a" }), @constCast(&[8]Unit{ si.m, si.kg, si.s, si.A, si.K, si.mol, si.cd, si.rad }));

pub const BASE_UNIT_CGS: SystemUnit = .init(@constCast(&[7][]const u8{ "l", "m", "t", "th", "n", "j", "a" }), @constCast(&[7]Unit{ si.cm, si.g, si.s, si.K, si.mol, si.cd, si.rad }));

pub const MAX_SYSTEM_LENGTH = 8;

pub const SystemUnit = struct {
    const Self = @This();

    len: usize = MAX_SYSTEM_LENGTH,
    dim_field: [MAX_SYSTEM_LENGTH][]const u8 = undefined,
    unit: [MAX_SYSTEM_LENGTH]Unit = undefined,

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

    fn getIdx(self: Self, field: []const u8) ?usize {
        if (field.len > 2) return null;
        for (0..self.len) |i| {
            if (std.mem.eql(u8, field, self.dim_field[i])) return i;
        }
        return null;
    }

    pub fn get(self: Self, field: []const u8) ?Unit {
        if (self.getIdx(field)) |idx| {
            return self.unit[idx];
        }
        return null;
    }
};

pub const System = enum {
    const Self = @This();

    SI,
    CGS,
    ESU,
    EMU,
    GAUSS,

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
