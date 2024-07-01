const std = @import("std");
const noise = @import("noise.zig");
const png = @import("png.zig");

fn usage() void {
	std.debug.print(
	\\Generates a png image of a cloud.
	\\Usage: cloud [options]
	\\    -s <seed>
	\\    -h <height>
	\\    -w <width>
	\\    -b <backround-color>
	\\    -f <foreground-color>
	\\    -o <outputfile>
	\\
	, .{});
}

pub fn main() !void {
    var seed: i64 = std.time.timestamp();
    var width: u64 = 512;
    var height: u64 = 512;
    var background: u24 = 0x00bbff;
    var foreground: u24 = 0xffffff;

    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const allocator = gpa.allocator();
	var filename = std.ArrayList(u8).init(allocator);
	defer filename.deinit();
	try filename.appendSlice("cloud.png");

    // Parse args without allocators
    // https://ziggit.dev/t/read-command-line-arguments/220/7
    var i: usize = 1;
    while (i < std.os.argv.len): (i += 1) {
        // seed
        if (std.mem.eql(u8, std.mem.span(std.os.argv[i]), "-s")) {
            if (std.os.argv.len > i) {
                seed = try std.fmt.parseInt(i64, std.mem.span(std.os.argv[i + 1]), 0);
                i += 1;
            } else {
            	usage();
            	return;
            }
        }
        // width
        else if (std.mem.eql(u8, std.mem.span(std.os.argv[i]), "-w")) {
            if (std.os.argv.len > i) {
                width = try std.fmt.parseInt(u64, std.mem.span(std.os.argv[i + 1]), 0);
                i += 1;
            } else {
            	usage();
            	return;
            }
        }
        // height
        else if (std.mem.eql(u8, std.mem.span(std.os.argv[i]), "-h")) {
            if (std.os.argv.len > i) {
                height = try std.fmt.parseInt(u64, std.mem.span(std.os.argv[i + 1]), 0);
                i += 1;
            } else {
            	usage();
            	return;
            }
        }
        // background color
        else if (std.mem.eql(u8, std.mem.span(std.os.argv[i]), "-b")) {
            if (std.os.argv.len > i) {
            	background = try std.fmt.parseInt(u24, std.mem.span(std.os.argv[i + 1]), 16);
            	i += 1;
            } else {
            	usage();
            	return;
            }
        }
        // foreground color
        else if (std.mem.eql(u8, std.mem.span(std.os.argv[i]), "-f")) {
            if (std.os.argv.len > i) {
            	foreground = try std.fmt.parseInt(u24, std.mem.span(std.os.argv[i + 1]), 16);
                i += 1;
            } else {
            	usage();
            	return;
            }
        }   
        // output file
        else if (std.mem.eql(u8, std.mem.span(std.os.argv[i]), "-o")) {
        	if (std.os.argv.len > i) {
        		filename.clearAndFree();
        		try filename.appendSlice(std.mem.span(std.os.argv[i + 1]));
        		i += 1;
        	} else {
        		usage();
        		return;
        	}
        } else {
            usage();
            return;
        }
    }
	std.debug.print("Creating {s} ({d}x{d}) with seed {d}, bg {x}, fg {x}.\n",
		.{filename.items, width, height, seed, background, foreground});

	var pixel: f32 = undefined;
	var image: png.pngStruct = try png.pngStruct.init(width, height, png.pngType.rgb, allocator);
    for (0..height) |h| {
    	for (0..width) |w| {
    		pixel = getNoise(seed, @as(f64, @floatFromInt(w)), @as(f64, @floatFromInt(h)));
			pixel -= @as(f32, @floatCast(getRadius(width, height, @as(u64, w), @as(u64, h))));
			const rgb = getColor(height, h, getRgb(background), getRgb(foreground), pixel);
			try image.setPixel(w, h, rgb);
		}
	}
	try image.build();
	try image.write(filename.items);
	try image.deinit();
}

fn getNoise(seed: i64, w: f64, h: f64) f32 {
	// feature size (24)
	const fsize = 48;
	// Four octaves: frequence N, N/2, N/4, N/8
	// with relatives amplitudes: 8:4:2:1
	const v0 = noise.noise2(seed, w / fsize / 8 , h / fsize / 8);
	const v1 = noise.noise2(seed, w / fsize / 4 , h / fsize / 4);
	const v2 = noise.noise2(seed, w / fsize / 2, h / fsize / 2);
	const v3 = noise.noise2(seed, w / fsize / 1, h / fsize / 1);
	const value = v0 * 8 / 13.0 + v1 * 4 / 13.0 + v2 * 2 / 13.0 + v3 * 1 / 13.0;
	// Three octaves: frequency N, N/2 and N/4
	// with relative amplitudes 4:2:1.
	//var value = v1 * 4 / 7.0 + v2 * 2 / 7.0 + v3 * 1 / 7.0;
	return (value + 1) / 2;
}

fn getRadius(width: u64, height: u64, w: u64, h: u64) f64 {
	const ch: f64 = @as(f64, @floatFromInt(height)) / 2;
	const cw: f64 = @as(f64, @floatFromInt(width)) / 2;
	const fw: f64 = @as(f64, @floatFromInt(w));
	const fh: f64 = @as(f64, @floatFromInt(h));
	var value = std.math.sqrt( (ch-fh)*(ch-fh) + (cw-fw)*(cw-fw) );
	value = value / (@as(f64, @floatFromInt(height)) / 1.5);
	if (value > 1) {
		value = 1;
	}
	return value;
}

fn getRgb(color: u24) [3]u8 {
	var rgb = [3]u8 {0, 0, 0};
	rgb[0] = @as(u8, @truncate(color >> 16 & 0xFF));
	rgb[1] = @as(u8, @truncate(color >> 8 & 0xFF));
	rgb[2] = @as(u8, @truncate(color & 0xFF));
	return rgb;
}

fn getColor(height: u64, h: u64, bg: [3]u8, fg: [3]u8, pixel: f32) u32 {
	var rgb = [4]u8 {0x00, 0x00, 0x00, 0x00};
	if (pixel <= 0) {
		rgb[0] = getGradient(height, h, bg[0]);
		rgb[1] = getGradient(height, h, bg[1]);
		rgb[2] = getGradient(height, h, bg[2]);
	} else {
		rgb[0] = getFgGradient(height, h, bg[0], fg[0], pixel);
		rgb[1] = getFgGradient(height, h, bg[1], fg[1], pixel);
		rgb[2] = getFgGradient(height, h, bg[2], fg[2], pixel);
	}
	if (rgb[0] > 255) {rgb[0] = 255;}
	if (rgb[0] < 0) {rgb[0] = 0;}
	if (rgb[1] > 255) {rgb[1] = 255;}
	if (rgb[1] < 0) {rgb[1] = 0;}
	if (rgb[2] > 255) {rgb[2] = 255;}
	if (rgb[2] < 0) {rgb[2] = 0;}

	return std.mem.readIntBig(u32, rgb[0..4]);
}

fn getGradient(height: u64, h: u64, color: u8) u8 {
	// we want a gradient background, darker on top.
	// This calculates the gradient color for a given row.
	return @as(u8, @intFromFloat(@floor(
		@as(f32, @floatFromInt(color)) +
		@as(f32, @floatFromInt(255-color)) * (@as(f64, @floatFromInt(h)) / @as(f64, @floatFromInt(height)) / 5.0)
	) ));
}

fn getFgGradient(height: u64, h: u64, bg: u8, fg: u8, pixel: f32) u8 {
	// default formula is r = floor( noise*r);
	// but we want our background color instead of black,
	// and to transition from background to foreground.
	// We could also transition from background to white:
	// r = floor (getGradient(height, h, br) * (255 - getGradient(height, h, fr) * 1.25);
	// https://stackoverflow.com/questions/21835739/smooth-color-transition-algorithm
	//rgb[0] = @floor( (1-pixel) * getGradient(height, h, bg[0]) + pixel * getGradient(height, h, fg[0]) );
	return @as(u8, @intFromFloat( @floor(
		(1 - pixel) * @as(f32, @floatFromInt(getGradient(height, h, bg) ))
		+ pixel * @as(f32, @floatFromInt(getGradient(height, h, fg) ))
	)));
}




test "simple test" {
    // do nothing
}
