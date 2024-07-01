### Cloud: generates clouds

Creates a random image of a cloud from a procedurally generated noise.

![A generated cloud](https://github.com/pjfichet/cloud/raw/master/cloud.png?raw=true)

### Credits and Licenses

The noise library (src/noise.zig) is a port of KdotJPG OpenSimplex2 noise.
- https://github.com/KdotJPG/OpenSimplex2
- public domain

The png library (src/png.zig) is a port of Michael Schwars libattopng.
- https://github.com/misc0110/libattopng
- Mit license

The main file (src/main.zig) is written from scratch.
- ISC license.

### Build and install

You will need the stable zig release (version 0.12).
To build and install in `$HOME/.local/bin/cloud`, run:

	zig build -Doptimize=ReleaseSafe -p $HOME/.local/

### usage

```
$ cloud [options]
    -s <seed> (default to std.time.timestamp())
    -w <width> (default to 512)
    -h <height> (default to 512)
    -b <backround-color> (default to 0088ee)
    -f <foreground-color> (default to ffffff)
    -o <outputfile> (default to cloud.png)
```
