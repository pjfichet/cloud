const std = @import("std");

const LIBATTOPNG_ADLER_BASE = 65521;

const crc32 = [_]u32 {
	    0x00000000, 0x77073096, 0xee0e612c, 0x990951ba, 0x076dc419, 0x706af48f, 0xe963a535, 0x9e6495a3, 0x0edb8832,
	    0x79dcb8a4, 0xe0d5e91e, 0x97d2d988, 0x09b64c2b, 0x7eb17cbd, 0xe7b82d07, 0x90bf1d91, 0x1db71064, 0x6ab020f2,
	    0xf3b97148, 0x84be41de, 0x1adad47d, 0x6ddde4eb, 0xf4d4b551, 0x83d385c7, 0x136c9856, 0x646ba8c0, 0xfd62f97a,
	    0x8a65c9ec, 0x14015c4f, 0x63066cd9, 0xfa0f3d63, 0x8d080df5, 0x3b6e20c8, 0x4c69105e, 0xd56041e4, 0xa2677172,
	    0x3c03e4d1, 0x4b04d447, 0xd20d85fd, 0xa50ab56b, 0x35b5a8fa, 0x42b2986c, 0xdbbbc9d6, 0xacbcf940, 0x32d86ce3,
	    0x45df5c75, 0xdcd60dcf, 0xabd13d59, 0x26d930ac, 0x51de003a, 0xc8d75180, 0xbfd06116, 0x21b4f4b5, 0x56b3c423,
	    0xcfba9599, 0xb8bda50f, 0x2802b89e, 0x5f058808, 0xc60cd9b2, 0xb10be924, 0x2f6f7c87, 0x58684c11, 0xc1611dab,
	    0xb6662d3d, 0x76dc4190, 0x01db7106, 0x98d220bc, 0xefd5102a, 0x71b18589, 0x06b6b51f, 0x9fbfe4a5, 0xe8b8d433,
	    0x7807c9a2, 0x0f00f934, 0x9609a88e, 0xe10e9818, 0x7f6a0dbb, 0x086d3d2d, 0x91646c97, 0xe6635c01, 0x6b6b51f4,
	    0x1c6c6162, 0x856530d8, 0xf262004e, 0x6c0695ed, 0x1b01a57b, 0x8208f4c1, 0xf50fc457, 0x65b0d9c6, 0x12b7e950,
	    0x8bbeb8ea, 0xfcb9887c, 0x62dd1ddf, 0x15da2d49, 0x8cd37cf3, 0xfbd44c65, 0x4db26158, 0x3ab551ce, 0xa3bc0074,
	    0xd4bb30e2, 0x4adfa541, 0x3dd895d7, 0xa4d1c46d, 0xd3d6f4fb, 0x4369e96a, 0x346ed9fc, 0xad678846, 0xda60b8d0,
	    0x44042d73, 0x33031de5, 0xaa0a4c5f, 0xdd0d7cc9, 0x5005713c, 0x270241aa, 0xbe0b1010, 0xc90c2086, 0x5768b525,
	    0x206f85b3, 0xb966d409, 0xce61e49f, 0x5edef90e, 0x29d9c998, 0xb0d09822, 0xc7d7a8b4, 0x59b33d17, 0x2eb40d81,
	    0xb7bd5c3b, 0xc0ba6cad, 0xedb88320, 0x9abfb3b6, 0x03b6e20c, 0x74b1d29a, 0xead54739, 0x9dd277af, 0x04db2615,
	    0x73dc1683, 0xe3630b12, 0x94643b84, 0x0d6d6a3e, 0x7a6a5aa8, 0xe40ecf0b, 0x9309ff9d, 0x0a00ae27, 0x7d079eb1,
	    0xf00f9344, 0x8708a3d2, 0x1e01f268, 0x6906c2fe, 0xf762575d, 0x806567cb, 0x196c3671, 0x6e6b06e7, 0xfed41b76,
	    0x89d32be0, 0x10da7a5a, 0x67dd4acc, 0xf9b9df6f, 0x8ebeeff9, 0x17b7be43, 0x60b08ed5, 0xd6d6a3e8, 0xa1d1937e,
	    0x38d8c2c4, 0x4fdff252, 0xd1bb67f1, 0xa6bc5767, 0x3fb506dd, 0x48b2364b, 0xd80d2bda, 0xaf0a1b4c, 0x36034af6,
	    0x41047a60, 0xdf60efc3, 0xa867df55, 0x316e8eef, 0x4669be79, 0xcb61b38c, 0xbc66831a, 0x256fd2a0, 0x5268e236,
	    0xcc0c7795, 0xbb0b4703, 0x220216b9, 0x5505262f, 0xc5ba3bbe, 0xb2bd0b28, 0x2bb45a92, 0x5cb36a04, 0xc2d7ffa7,
	    0xb5d0cf31, 0x2cd99e8b, 0x5bdeae1d, 0x9b64c2b0, 0xec63f226, 0x756aa39c, 0x026d930a, 0x9c0906a9, 0xeb0e363f,
	    0x72076785, 0x05005713, 0x95bf4a82, 0xe2b87a14, 0x7bb12bae, 0x0cb61b38, 0x92d28e9b, 0xe5d5be0d, 0x7cdcefb7,
	    0x0bdbdf21, 0x86d3d2d4, 0xf1d4e242, 0x68ddb3f8, 0x1fda836e, 0x81be16cd, 0xf6b9265b, 0x6fb077e1, 0x18b74777,
	    0x88085ae6, 0xff0f6a70, 0x66063bca, 0x11010b5c, 0x8f659eff, 0xf862ae69, 0x616bffd3, 0x166ccf45, 0xa00ae278,
	    0xd70dd2ee, 0x4e048354, 0x3903b3c2, 0xa7672661, 0xd06016f7, 0x4969474d, 0x3e6e77db, 0xaed16a4a, 0xd9d65adc,
	    0x40df0b66, 0x37d83bf0, 0xa9bcae53, 0xdebb9ec5, 0x47b2cf7f, 0x30b5ffe9, 0xbdbdf21c, 0xcabac28a, 0x53b39330,
	    0x24b4a3a6, 0xbad03605, 0xcdd70693, 0x54de5729, 0x23d967bf, 0xb3667a2e, 0xc4614ab8, 0x5d681b02, 0x2a6f2b94,
	    0xb40bbe37, 0xc30c8ea1, 0x5a05df1b, 0x2d02ef8d
};

fn swap32(num: u32) u32 {
	return ((num >> 24) & 0xff) |
	       ((num << 8) & 0xff0000) |
	       ((num >> 8) & 0xff00) |
	       ((num << 24) & 0xff000000);
}

fn crc(data: []const u8, len: usize, oldcrc: u32) u32 {
	var newcrc = oldcrc;
	for (0..len) |i| {
	    newcrc = crc32[(newcrc ^ data[i]) & 255] ^ (newcrc >> 8);
		//std.debug.print("data: {d}, newcrc: {d}\n", .{data[i], newcrc});
	}
	//std.debug.print("data: {s}, crc: {d}\n", .{data, newcrc});
	return newcrc;
}

pub const pngType = enum(u8) {
	grayscale = 0,          // 256 shades of gray, 8bit per pixel
	rgb = 2,                // 24bit RGB values
	palette = 3,            // Up to 256 RGBA palette colors, 8bit per pixel
	grayscale_alpha = 4,    // 256 shades of gray plus alpha channel, 16bit per pixel
	rgba = 6                // 24bit RGB values plus 8bit alpha channel
};

pub const pngStruct = struct {
	allocator: std.mem.Allocator = undefined,
	pngtype: pngType = .rgb,         // File type
	capacity: usize = 0,             // Reserved memory for raw data
	data: []u8 = undefined,          // Raw pixel data, format depends on type
	palette: []u32 = undefined,      // Palette for image
	palette_length: usize = 0,       // Entries for palette, 0 if unused
	width: usize = 0,                // Image width
	height: usize = 0,               // Image height

	out: []u8 = undefined,           // Buffer to store final PNG
	out_pos: usize = 0,              // Current size of output buffer
	out_capacity: usize = 0,         // Capacity of output buffer
	crc: u32 = 0,                    // Currecnt CRC32 checksum
	s1: u16 = 1,                     // Helper variables for Adler checksum
	s2: u16 = 0,                     // Helper variables for Adler checksum
	bpp: usize = 0,                  // Bytes per pixel


	pub fn init(width: usize, height: usize, pngtype: pngType, allocator: std.mem.Allocator) !pngStruct {
		var png = pngStruct {
			.allocator = allocator,
			.width = width,
	    	.height = height,
	    	.capacity = width * height,
		    .pngtype = pngtype,
		};
		switch (png.pngtype) {
	    	.grayscale => png.bpp = 1,
	    	.rgb => {png.bpp = 3; png.capacity *= 4;},
	    	.palette => png.bpp = 1,
	    	.grayscale_alpha => {png.bpp = 2; png.capacity *=2 ;},
	    	.rgba => {png.bpp = 4; png.capacity *= 4;},
	    }
		png.data = try png.allocator.alloc(u8, png.capacity);
		errdefer png.allocator.free(png.data);
		if (png.pngtype == .palette) {
			png.palette = try png.allocator.alloc(u32, 256);
			errdefer png.allocator.free(png.palette);
		}
	    return png;
	}

	pub fn deinit(png: *pngStruct) !void {
		if (png.pngtype == .palette) {
			png.allocator.free(png.palette);
		}
		if (png.out_capacity > 0) {
			png.allocator.free(png.out);
		}
		png.allocator.free(png.data);
	}

	pub fn setPalette(png: *pngStruct, palette: []const u32) !void {
		if (palette.len > 256) {
			return error.PaletteTooBig;
		}
		@memcpy(png.palette, palette);
		png.palette_length = palette.len;
	}

	pub fn setPixel(png: *pngStruct, x: usize, y: usize, color: u32) !void {
	    if (x >= png.width or y >= png.height) {
	    	return error.WrongCoordinates;
	    }
	    // https://stackoverflow.com/questions/74701099/how-to-insert-a-u32-into-a-u8-array-in-zig
	    if (png.pngtype == .palette or png.pngtype == .grayscale) {
	        png.data[x + y * png.width] = @as(u8, @intCast(color & 0xff));
	    } else if (png.pngtype == .grayscale_alpha) {
	    	// writeIntBig, writeIntNative
	    	std.mem.writeIntBig(u16, png.data[(x + y * png.width) * @sizeOf(u16)..][0..@sizeOf(u16)], @as(u16, @intCast(color & 0xffff)) );
	    } else {
	    	//std.debug.print("index({}, {}, {}): {}\n", .{x, y, png.width, x + y * png.width * @sizeOf(u32)});
	    	std.mem.writeIntBig(u32, png.data[(x + y * png.width) * @sizeOf(u32)..][0..@sizeOf(u32)], color);
	    }
	}

	pub fn getPixel(png: *pngStruct, x: usize, y: usize) u32 {
	    var pixel: u32 = 0;
	    if (x >= png.width or y >= png.height) {
	        return pixel;
	    }
	    if (png.pngtype == .palette or png.pngtype == .grayscale) {
	        //pixel = (uint32_t)(png.data[x + y * png.width] & 0xff);
	        pixel = @as(u32, png.data[x + y * png.width] & 0xff);
	    } else if (png.pngtype == .grayscale_alpha) {
	        //pixel = (uint32_t)(((uint16_t *) png.data)[x + y * png.width] & 0xffff);
	        pixel = @as(u32, std.mem.readIntBig(u16, png.data[(x + y * png.width) * @sizeOf(u16)..][0..@sizeOf(u16)]) & 0xffff);
	    } else {
	        //pixel = ((uint32_t *) png.data)[x + y * png.width];
	        pixel = std.mem.readIntBig(u32, png.data[(x + y * png.width) * @sizeOf(u32)..][0..@sizeOf(u32)]);
	    }
	    return pixel;
	}

	/////////////////////////////////////////////////////////////////////////
	// Utils

	fn outRawWrite(png: *pngStruct, data: []const u8, len: usize) void {
	    for (0..len) |i| {
	        png.out[png.out_pos] = data[i];
	        png.out_pos += 1;
		}
	}

	fn outRawU32(png: *pngStruct, val: u32) void {
		std.mem.writeIntNative(u32, png.out[png.out_pos..][0..@sizeOf(u32)], val);
	    png.out_pos += 4;
	}

	fn outRawU16(png: *pngStruct, val: u16) void {
	    std.mem.writeIntNative(u16, png.out[png.out_pos..][0..@sizeOf(u16)], val);
	    png.out_pos += 2;
	}

	fn outRawU8(png: *pngStruct, val: u8) void {
		png.out[png.out_pos] = val;
	    png.out_pos += 1;
	}

	fn newChunk(png: *pngStruct,name: []const u8, len: usize) void {
	    png.crc = 0xffffffff;
	    png.outRawU32(swap32(@as(u32, @intCast(len))) );
	    png.crc = crc(name, 4, png.crc);
	    png.outRawWrite(name, 4);
	}

	fn endChunk(png: *pngStruct) void {
	    png.outRawU32(swap32(~png.crc));
	    //std.debug.print("end chunk crc: {d} {d} {d}\n", .{png.crc, ~png.crc, swap32(~png.crc)});
	}

	fn outU32(png: *pngStruct, val: u32) void {
		var data: [4]u8 = undefined;
		std.mem.writeIntNative(u32, data[0..4], val);
	    png.crc = crc(&data, 4, png.crc);
	    png.outRawU32(val);
	}

	fn outU16(png: *pngStruct, val: u16) void {
		var data: [2]u8 = undefined;
		std.mem.writeIntNative(u16, data[0..2], val);
	    png.crc = crc(&data, 2, png.crc);
	    png.outRawU16(val);
	}

	fn outU8(png: *pngStruct, val: u8) void {
		const data = [_]u8 {val};
	    png.crc = crc(&data, 1, png.crc);
	    png.outRawU8(val);
	}

	fn outWrite(png: *pngStruct, data: []const u8, len: usize) void {
		//todo: png->crc = libattopng_crc((const unsigned char *) data, len, png->crc);
	    png.crc = crc(data, len, png.crc);
	    png.outRawWrite(data, len);
	}

	fn outWriteAdler(png: *pngStruct, val: u8) void {
	    //png.outWrite(png, (char *) &data, 1);
		const data = [_]u8 {val};
	    png.outWrite(&data, 1);
		//std.debug.print("s1: {}, s2: {}\n", .{png.s1, png.s2});
	    //png.s1 = @as(u16, (png.s1 +% val) % LIBATTOPNG_ADLER_BASE);
	    //png.s2 = @as(u16, (png.s2 + png.s1) % LIBATTOPNG_ADLER_BASE);
		// modulo should ensure we don't need to wrap. But we need to handle
		// integer overflow.
		const s1: u32 = @as(u32, png.s1) + @as(u32, val);
		png.s1 = @as(u16, @intCast(s1 % LIBATTOPNG_ADLER_BASE));
		const s2: u32 = @as(u32, png.s2)  + @as(u32, png.s1);
		png.s2 = @as(u16, @intCast(s2 % LIBATTOPNG_ADLER_BASE));

	    //png.s1 = (png.s1 +% val) % LIBATTOPNG_ADLER_BASE;
	    //png.s2 = (png.s2 +% png.s1) % LIBATTOPNG_ADLER_BASE;
		//std.debug.print("s1: {}, s2: {}\n", .{png.s1, png.s2});
	}

	fn pixelHeader(png: *pngStruct, offset: usize, bpl: usize) void {
	    if (offset > bpl) {
	        // not the last line
	        //libattopng_outWrite(png, "\0", 1);
	        const nul = [1]u8 {0};
	        png.outWrite(&nul, 1);
	        png.outU16(@as(u16, @intCast(bpl)) );
	        png.outU16(~@as(u16, @intCast(bpl)) );
	    } else {
	        // last line
	        //libattopng_outWrite(png, "\1", 1);
	        const one = [1]u8 {1};
	        png.outWrite(&one, 1);
	        png.outU16(@as(u16, @intCast(offset)) );
	        png.outU16(~@as(u16, @intCast(offset)) );
	    }
	}

	////////////////////////////////////////////////////////////////////////////
	pub fn build(png: *pngStruct) !void {
	    // size_t index, bpl, raw_size, size, p, pos, corr;
	    //unsigned char *pixel;
	    if (png.out_capacity > 0) {
	        // delete old output if any
	        png.allocator.free(png.out);
	    }
	    png.out_capacity = png.capacity + 4096 * 8 + png.width * png.height;
	    png.out = try png.allocator.alloc(u8, png.out_capacity);
	    errdefer png.allocator.free(png.out);

		// PNG header
	    // libattopng_outRawWrite(png, "\211PNG\r\n\032\n", 8);
		png.outRawU8(137);
		png.outRawWrite("PNG\r\n", 5);
		png.outRawU8(26);
		png.outRawWrite("\n", 1);

		// IHDR
		png.newChunk("IHDR", 13);
		png.outU32(swap32( @as(u32, @intCast(png.width)) ));
		png.outU32(swap32(@as(u32, @intCast(png.height)) ));
	    png.outU8(8); // bit depth
	    png.outU8(@intFromEnum(png.pngtype));
	    png.outU8(0); // compression
	    png.outU8(0); // filter
	    png.outU8(0); // interlace method
	    png.endChunk();


	    // palette //
	    if (png.pngtype == .palette) {
	        var entry: [3]u8 = undefined;
	        var s: usize = png.palette_length;
	        if (s < 16) {
	            s = 16; // minimum palette length
	        }
	        png.newChunk("PLTE", 3 * s);
	        for (0..s) |i| {
	            entry[0] = @as(u8, @intCast(png.palette[i] & 255));
	            entry[1] = @as(u8, @intCast((png.palette[i] >> 8) & 255));
	            entry[2] = @as(u8, @intCast((png.palette[i] >> 16) & 255));
	            png.outWrite(&entry, 3);
	        }
	        png.endChunk();

	        // transparency
	        png.newChunk("tRNS", s);
	        for (0..s) |j| {
	            entry[0] = @as(u8, @intCast((png.palette[j] >> 24) & 255));
	            png.outWrite(&entry, 1);
	        }
	        png.endChunk();
	    }

	    // data
	    var bpl = 1 + png.bpp * png.width;
	    if (bpl >= 65536) {
	        std.debug.print("[libattopng] ERROR: maximum supported width for this type of PNG is {d} pixel\n",
	        	.{65535 / @as(u32, @intCast(png.bpp))});
	        return error.pngTooBig;
	    }
	    var raw_size = png.height * bpl;
	    var size = 2 + png.height * (5 + bpl) + 4;
	    //std.debug.print("width: {d}, height: {d}, bpp: {d}, bpl: {d}\n", .{png.width, png.height, png.bpp, bpl});
	    png.newChunk("IDAT", size);
	    //png.outWrite(png, "\170\332", 2);
	    png.outU8(120);
	    png.outU8(218);

	    var pixel: usize = 0;
	    var index: usize = 0;
	    var pos: usize = 0;
	    while (pos < png.width * png.height) : (pos += 1) {
	        if (index == 0) {
	            // line header
	            png.pixelHeader(raw_size, bpl);
	            png.outWriteAdler(0); // no filter
	            raw_size -= 1;
	        }

	        // pixel
	        for (0..png.bpp) |_| {
	        	//std.debug.print("pixel[{d}]: {x}\n", .{pixel, png.data[pixel]});
	            png.outWriteAdler(png.data[pixel]);
	            pixel += 1;
	        }
	        if (png.pngtype == .rgb) {
				pixel += 1;
			}

	        raw_size -= png.bpp;
	        index = (index + 1) % png.width;
	    }
	    // checksum
	    //std.debug.print("s1: {}, s2: {}\n", .{png.s1, png.s2});
	    png.s1 %= LIBATTOPNG_ADLER_BASE;
	    png.s2 %= LIBATTOPNG_ADLER_BASE;
	    //std.debug.print("s1: {}, s2: {}\n", .{png.s1, png.s2});
	    png.outU32(swap32(@as(u32, (@as(u32, png.s2) << 16) | png.s1)) );
	    png.endChunk();

	    // end of image
	    png.newChunk("IEND", 0);
	    png.endChunk();
	}

	pub fn write(png: *pngStruct, filename: []const u8) !void {
		const file = try std.fs.cwd().createFile(filename, .{});
		defer file.close();
		const bytes_written = try file.write(png.out[0..png.out_pos]);
		_ = bytes_written;
	}

}; // end struct 

////////////////////////////////////////////////////////////////////

test "pngstruct" {
	std.debug.print("Start...\n", .{});
	const colors = [_]u32 {
		0x00000000, 0x11111111, 0x22222222, 0x33333333,
		0x44444444, 0x55555555, 0x66666666, 0x77777777,
		0x88888888, 0x99999999, 0xaaaaaaaa, 0xbbbbbbbb,
		0xcccccccc, 0xdddddddd, 0xeeeeeeee, 0xffffffff
	};
	var png: pngStruct = try pngStruct.init(4, 4, .rgb, std.testing.allocator);
	var pixel: u32 = undefined;
	for (0..4) |x| {
		for (0..4) |y| {
			try png.setPixel(x, y, colors[x  + y * 4]);
			pixel = png.getPixel(x, y);
//			std.debug.print("pixel {d} {d}: {x}\n", .{x, y, pixel});
			try std.testing.expectEqual(pixel, colors[x + y * 4]);
		}
	}
	//std.debug.print("data: {any}\n", .{png.data});
//	for (0..png.data.len) |z| {
//		std.debug.print("data[{d}] = {x}\n", .{z, png.data[z]});
//	}
//	for (0..16) |p| {
//		pixel = std.mem.readIntNative(u32, png.data[p..][0..@sizeOf(u32)]);
//		std.debug.print("pixel[{d}] = {x}\n", .{p, pixel});
//	}
	try png.build();
	try png.write("test.png");
	try png.deinit();
}

test "color" {
	var png: pngStruct = try pngStruct.init(4, 4, .rgb, std.testing.allocator);
	for (0..4) |x| {
		for (0..4) |y| {
			try png.setPixel(x, y, 0x00bbff00);
		}
	}
	try png.build();
	try png.write("test.png");
	try png.deinit();
}

test "swap32" {
	const x = 0xaabbccdd;
	std.debug.print("swap32({x}) = {x}\n", .{x, swap32(x)});
}
