const std = @import("std");

pub fn expectApproxEqAbsIter(expected: anytype, actual: anytype, tolerance: anytype) !void {
    try std.testing.expectEqual(expected.len, actual.len);
    for (expected, actual) |e, a| try std.testing.expectApproxEqAbs(e, a, tolerance);
}

pub fn expectApproxEqRelIter(expected: anytype, actual: anytype, tolerance: anytype) !void {
    try std.testing.expectEqual(expected.len, actual.len);
    for (expected, actual) |e, a| try std.testing.expectApproxEqRel(e, a, tolerance);
}
