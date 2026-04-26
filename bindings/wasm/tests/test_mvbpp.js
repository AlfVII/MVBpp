/**
 * mvbpp WASM binding tests.
 *
 * Run with:  node bindings/wasm/tests/test_mvbpp.js
 * Requires:  bindings/wasm/dist/mvbpp.js  (built via `npm run build`)
 */

'use strict';

const fs   = require('fs');
const path = require('path');

const DIST_JS   = path.resolve(__dirname, '../dist/mvbpp.js');
const TESTS_DIR = __dirname;

// ── tiny test harness ─────────────────────────────────────────────────────────

let passed = 0, failed = 0;
const results = [];

function test(name, fn) {
    try {
        fn();
        console.log(`  ✓ ${name}`);
        passed++;
    } catch (e) {
        console.error(`  ✗ ${name}`);
        console.error(`    ${e.message}`);
        failed++;
    }
}

function assert(cond, msg) {
    if (!cond) throw new Error(msg || 'Assertion failed');
}

function assertEqual(a, b, msg) {
    if (a !== b) throw new Error(msg || `Expected ${b}, got ${a}`);
}

// ── helpers ───────────────────────────────────────────────────────────────────

function loadMagneticJson(filename) {
    const raw  = fs.readFileSync(path.join(TESTS_DIR, filename), 'utf8');
    const data = JSON.parse(raw);
    return JSON.stringify(data.magnetic);
}

function isStepBytes(buf) {
    // STEP files begin with "ISO-10303-21"
    const header = Buffer.from(buf.slice(0, 12));
    return header.toString('ascii') === 'ISO-10303-21';
}

function bufferToString(buf) {
    return Buffer.from(buf).toString('latin1');
}

// ── main ──────────────────────────────────────────────────────────────────────

async function main() {
    if (!fs.existsSync(DIST_JS)) {
        console.error(`ERROR: ${DIST_JS} not found. Run 'npm run build' in bindings/wasm/ first.`);
        process.exit(1);
    }

    const createMvbpp = require(DIST_JS);
    const mvbpp = await createMvbpp();

    console.log('\nmvbpp WASM — drawMagneticToBuffer');
    console.log('─'.repeat(50));

    // ── drawMagneticToBuffer ─────────────────────────────────────────────────

    test('returns Uint8Array for concentric_basic', () => {
        const j      = loadMagneticJson('concentric_basic.json');
        const result = mvbpp.drawMagneticToBuffer(j);
        assert(result instanceof Uint8Array, 'Expected Uint8Array');
    });

    test('output is non-empty', () => {
        const j      = loadMagneticJson('concentric_basic.json');
        const result = mvbpp.drawMagneticToBuffer(j);
        assert(result.length > 1000, `Too small: ${result.length} bytes`);
    });

    test('output has valid STEP header', () => {
        const j      = loadMagneticJson('concentric_basic.json');
        const result = mvbpp.drawMagneticToBuffer(j);
        assert(isStepBytes(result), 'Does not start with ISO-10303-21');
    });

    test('STEP contains solid geometry markers', () => {
        const j    = loadMagneticJson('concentric_basic.json');
        const text = bufferToString(mvbpp.drawMagneticToBuffer(j));
        assert(
            text.includes('CLOSED_SHELL') || text.includes('ADVANCED_BREP'),
            'No solid geometry markers found in STEP output'
        );
    });

    test('ETD49 produces valid STEP', () => {
        const j      = loadMagneticJson('ETD49_N87_10uH_5T.json');
        const result = mvbpp.drawMagneticToBuffer(j);
        assert(isStepBytes(result), 'ETD49 STEP header invalid');
        assert(result.length > 10_000, `ETD49 output too small: ${result.length}`);
    });

    test('ETD49 produces more geometry than concentric_basic', () => {
        const basic = mvbpp.drawMagneticToBuffer(loadMagneticJson('concentric_basic.json'));
        const etd   = mvbpp.drawMagneticToBuffer(loadMagneticJson('ETD49_N87_10uH_5T.json'));
        assert(etd.length > basic.length,
            `ETD49 (${etd.length}) should be larger than basic (${basic.length})`);
    });

    // ── drawMagnetic (to virtual FS path) ────────────────────────────────────

    console.log('\nmvbpp WASM — drawMagnetic (virtual FS)');
    console.log('─'.repeat(50));

    test('writes file to emscripten virtual FS', () => {
        const j    = loadMagneticJson('concentric_basic.json');
        const out  = '/tmp/test_basic.step';
        mvbpp.drawMagnetic(j, out);
        const data = mvbpp.FS.readFile(out);
        assert(data.length > 1000, `Virtual FS file too small: ${data.length}`);
        mvbpp.FS.unlink(out);
    });

    test('virtual FS file has valid STEP header', () => {
        const j   = loadMagneticJson('concentric_basic.json');
        const out = '/tmp/test_header.step';
        mvbpp.drawMagnetic(j, out);
        const data = mvbpp.FS.readFile(out);
        assert(isStepBytes(data), 'Virtual FS file does not start with ISO-10303-21');
        mvbpp.FS.unlink(out);
    });

    test('virtual FS file matches buffer output', () => {
        const j   = loadMagneticJson('concentric_basic.json');
        const out = '/tmp/test_match.step';
        mvbpp.drawMagnetic(j, out);
        const fileData   = mvbpp.FS.readFile(out);
        const bufferData = mvbpp.drawMagneticToBuffer(j);
        mvbpp.FS.unlink(out);
        const ratio = Math.abs(fileData.length - bufferData.length) /
                      Math.max(fileData.length, bufferData.length);
        assert(ratio < 0.01, `File vs buffer size differs by ${(ratio * 100).toFixed(1)}%`);
    });

    // ── DrawConfig ────────────────────────────────────────────────────────────

    console.log('\nmvbpp WASM — DrawConfig');
    console.log('─'.repeat(50));

    test('DrawConfig default produces valid STEP', () => {
        const j   = loadMagneticJson('concentric_basic.json');
        const cfg = new mvbpp.DrawConfig();
        const res = mvbpp.drawMagneticToBufferWithConfig(j, cfg);
        assert(isStepBytes(res), 'DrawConfig default: invalid STEP header');
        cfg.delete();
    });

    test('no bobbin produces smaller output', () => {
        const j = loadMagneticJson('concentric_basic.json');

        const cfgWith = new mvbpp.DrawConfig();
        cfgWith.includeBobbin = true;
        const withBobbin = mvbpp.drawMagneticToBufferWithConfig(j, cfgWith);
        cfgWith.delete();

        const cfgWithout = new mvbpp.DrawConfig();
        cfgWithout.includeBobbin = false;
        const withoutBobbin = mvbpp.drawMagneticToBufferWithConfig(j, cfgWithout);
        cfgWithout.delete();

        assert(withBobbin.length > withoutBobbin.length,
            `With bobbin (${withBobbin.length}) should be larger than without (${withoutBobbin.length})`);
    });

    test('scale=1000 accepted without error', () => {
        const j   = loadMagneticJson('concentric_basic.json');
        const cfg = new mvbpp.DrawConfig();
        cfg.scale = 1000.0;
        const res = mvbpp.drawMagneticToBufferWithConfig(j, cfg);
        assert(isStepBytes(res), 'scale=1000 produced invalid STEP');
        cfg.delete();
    });

    // ── error handling ────────────────────────────────────────────────────────

    console.log('\nmvbpp WASM — error handling');
    console.log('─'.repeat(50));

    test('invalid JSON throws', () => {
        let threw = false;
        try { mvbpp.drawMagneticToBuffer('not json'); } catch { threw = true; }
        assert(threw, 'Expected exception for invalid JSON');
    });

    test('empty object throws', () => {
        let threw = false;
        try { mvbpp.drawMagneticToBuffer('{}'); } catch { threw = true; }
        assert(threw, 'Expected exception for empty object');
    });

    // ── summary ───────────────────────────────────────────────────────────────

    console.log('\n' + '─'.repeat(50));
    console.log(`${passed + failed} tests: ${passed} passed, ${failed} failed`);

    if (failed > 0) process.exit(1);
}

main().catch(e => { console.error(e); process.exit(1); });
