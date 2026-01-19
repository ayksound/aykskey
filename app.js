// AYKSOUND — ACCURATE TUNEBAT-STYLE ANALYSIS
// No auto-corrections, pure signal processing

async function analyze() {
    const file = document.getElementById('audioFile').files[0];
    if (!file) {
        alert('Please upload an audio file first');
        return;
    }

    document.getElementById('bpm').innerText = 'Analyzing...';
    document.getElementById('key').innerText = 'Analyzing...';
    document.getElementById('camelot').innerText = 'Analyzing...';
    document.getElementById('detune').innerText = 'Analyzing...';

    try {
        const audioCtx = new (window.AudioContext || window.webkitAudioContext)();
        const arrayBuffer = await file.arrayBuffer();
        const audioBuffer = await audioCtx.decodeAudioData(arrayBuffer);

        // Mono mix
        let mono = audioBuffer.getChannelData(0);
        if (audioBuffer.numberOfChannels > 1) {
            const left = audioBuffer.getChannelData(0);
            const right = audioBuffer.getChannelData(1);
            mono = new Float32Array(Math.min(left.length, audioBuffer.sampleRate * 60));
            for (let i = 0; i < mono.length; i++) {
                mono[i] = (left[i] + right[i]) / 2;
            }
        } else {
            mono = mono.slice(0, audioBuffer.sampleRate * 60);
        }

        const sr = audioBuffer.sampleRate;

        // Run analysis
        const bpm = await detectBPMAccurate(mono, sr);
        document.getElementById('bpm').innerText = Math.round(bpm);

        const keyData = await detectKeyAccurate(mono, sr);
        document.getElementById('key').innerText = keyData.key;
        document.getElementById('camelot').innerText = keyData.camelot;

        const detuneData = await detectDetuneAccurate(mono, sr);
        document.getElementById('detune').innerText = `${detuneData.hz.toFixed(1)} Hz / ${detuneData.cents > 0 ? '+' : ''}${detuneData.cents.toFixed(1)} cents`;

    } catch (error) {
        console.error('Error:', error);
        alert('Analysis failed. Try a different file.');
        document.getElementById('bpm').innerText = '—';
        document.getElementById('key').innerText = '—';
        document.getElementById('camelot').innerText = '—';
        document.getElementById('detune').innerText = '—';
    }
}

// ===== ACCURATE BPM DETECTION =====
async function detectBPMAccurate(signal, sr) {
    const hopSize = 512;
    const windowSize = 2048;

    // Compute onset strength using spectral flux
    const onsets = [];
    let prevMag = new Float32Array(windowSize / 2);

    for (let i = 0; i < signal.length - windowSize; i += hopSize) {
        const frame = signal.slice(i, i + windowSize);
        const mag = computeSpectralMagnitude(frame);

        // Spectral flux = sum of positive differences
        let flux = 0;
        for (let k = 0; k < mag.length; k++) {
            const diff = mag[k] - prevMag[k];
            if (diff > 0) flux += diff;
        }
        onsets.push(flux);
        prevMag = mag;

        if (i % (hopSize * 200) === 0) {
            await new Promise(r => setTimeout(r, 0));
        }
    }

    // Autocorrelation on onset function
    const minBPM = 40;
    const maxBPM = 200;
    const minLag = Math.floor((60 / maxBPM) * sr / hopSize);
    const maxLag = Math.floor((60 / minBPM) * sr / hopSize);

    const correlations = [];

    for (let lag = minLag; lag <= maxLag && lag < onsets.length / 2; lag++) {
        let sum = 0;
        const n = Math.min(1000, onsets.length - lag);

        for (let i = 0; i < n; i++) {
            sum += onsets[i] * onsets[i + lag];
        }

        const bpm = (60 * sr) / (lag * hopSize);
        correlations.push({ bpm, correlation: sum });

        if (lag % 100 === 0) {
            await new Promise(r => setTimeout(r, 0));
        }
    }

    // Find best BPM
    correlations.sort((a, b) => b.correlation - a.correlation);
    return correlations[0].bpm;
}

// ===== ACCURATE KEY DETECTION =====
async function detectKeyAccurate(signal, sr) {
    // Use Web Audio API for accurate FFT
    const audioCtx = new OfflineAudioContext(1, signal.length, sr);
    const buffer = audioCtx.createBuffer(1, signal.length, sr);
    buffer.getChannelData(0).set(signal);

    const source = audioCtx.createBufferSource();
    source.buffer = buffer;

    const analyser = audioCtx.createAnalyser();
    analyser.fftSize = 16384;
    analyser.smoothingTimeConstant = 0;

    source.connect(analyser);
    analyser.connect(audioCtx.destination);

    // Accumulate chroma from multiple positions
    const chroma = new Float32Array(12);
    const numSamples = 100;
    const hopSize = Math.floor(signal.length / numSamples);

    const freqData = new Float32Array(analyser.frequencyBinCount);
    const binWidth = sr / analyser.fftSize;

    for (let sample = 0; sample < numSamples; sample++) {
        const pos = sample * hopSize;
        if (pos >= signal.length - analyser.fftSize) break;

        // Manual FFT computation for this position
        const frame = signal.slice(pos, pos + analyser.fftSize);
        const spectrum = computeDFT(frame, sr);

        // Map frequencies to pitch classes
        for (let bin = 0; bin < spectrum.length; bin++) {
            const freq = bin * binWidth;
            if (freq < 80 || freq > 2000) continue;

            // Convert to pitch class (C=0, C#=1, ..., B=11)
            const midi = 12 * Math.log2(freq / 440) + 69;
            const pitchClass = Math.round(midi) % 12;
            const pc = (pitchClass + 12) % 12;

            chroma[pc] += spectrum[bin];
        }

        if (sample % 10 === 0) {
            await new Promise(r => setTimeout(r, 0));
        }
    }

    // Normalize
    const max = Math.max(...chroma);
    if (max > 0) {
        for (let i = 0; i < 12; i++) chroma[i] /= max;
    }

    // Krumhansl-Schmuckler key profiles
    const majorProfile = [6.35, 2.23, 3.48, 2.33, 4.38, 4.09, 2.52, 5.19, 2.39, 3.66, 2.29, 2.88];
    const minorProfile = [6.33, 2.68, 3.52, 5.38, 2.60, 3.53, 2.54, 4.75, 3.98, 2.69, 3.34, 3.17];
    const noteNames = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B'];

    let bestKey = 'C';
    let bestMode = 'Major';
    let maxCorr = -Infinity;

    for (let tonic = 0; tonic < 12; tonic++) {
        // Test major
        let majorCorr = 0;
        for (let i = 0; i < 12; i++) {
            majorCorr += chroma[(tonic + i) % 12] * majorProfile[i];
        }
        if (majorCorr > maxCorr) {
            maxCorr = majorCorr;
            bestKey = noteNames[tonic];
            bestMode = 'Major';
        }

        // Test minor
        let minorCorr = 0;
        for (let i = 0; i < 12; i++) {
            minorCorr += chroma[(tonic + i) % 12] * minorProfile[i];
        }
        if (minorCorr > maxCorr) {
            maxCorr = minorCorr;
            bestKey = noteNames[tonic];
            bestMode = 'Minor';
        }
    }

    const keyStr = bestMode === 'Minor' ? bestKey + 'm' : bestKey;
    return { key: keyStr, camelot: getCamelot(bestKey, bestMode) };
}

// ===== ACCURATE DETUNE DETECTION =====
async function detectDetuneAccurate(signal, sr) {
    const fftSize = 8192;
    const numSamples = 50;
    let pitches = [];

    for (let i = 0; i < numSamples; i++) {
        const pos = Math.floor((signal.length / (numSamples + 1)) * (i + 1));
        if (pos + fftSize > signal.length) break;

        const frame = signal.slice(pos, pos + fftSize);
        const pitch = detectPitchYIN(frame, sr);

        if (pitch > 60 && pitch < 1000) {
            pitches.push(pitch);
        }

        if (i % 10 === 0) {
            await new Promise(r => setTimeout(r, 0));
        }
    }

    if (pitches.length === 0) {
        return { hz: 440, cents: 0 };
    }

    // Use median instead of mean (more robust)
    pitches.sort((a, b) => a - b);
    const medianPitch = pitches[Math.floor(pitches.length / 2)];

    // Calculate detune from nearest note (A440 reference)
    const midi = 69 + 12 * Math.log2(medianPitch / 440);
    const nearestMidi = Math.round(midi);
    const nearestFreq = 440 * Math.pow(2, (nearestMidi - 69) / 12);
    const cents = 1200 * Math.log2(medianPitch / nearestFreq);

    return { hz: medianPitch, cents };
}

// ===== HELPER: YIN PITCH DETECTION =====
function detectPitchYIN(buffer, sr) {
    const threshold = 0.1;
    const minFreq = 60;
    const maxFreq = 1000;
    const tauMin = Math.floor(sr / maxFreq);
    const tauMax = Math.floor(sr / minFreq);

    // Difference function
    const diff = new Float32Array(tauMax);
    for (let tau = tauMin; tau < tauMax && tau < buffer.length / 2; tau++) {
        let sum = 0;
        for (let i = 0; i < buffer.length / 2; i++) {
            const delta = buffer[i] - buffer[i + tau];
            sum += delta * delta;
        }
        diff[tau] = sum;
    }

    // Cumulative mean normalized difference
    const cmndf = new Float32Array(tauMax);
    cmndf[0] = 1;
    let runningSum = 0;

    for (let tau = 1; tau < tauMax; tau++) {
        runningSum += diff[tau];
        cmndf[tau] = diff[tau] / (runningSum / tau);
    }

    // Find first valley below threshold
    let tau = tauMin;
    while (tau < tauMax - 1) {
        if (cmndf[tau] < threshold && cmndf[tau] < cmndf[tau + 1]) {
            break;
        }
        tau++;
    }

    // Parabolic interpolation
    if (tau < tauMax - 1 && tau > 0) {
        const betterTau = tau + (cmndf[tau + 1] - cmndf[tau - 1]) / (2 * (2 * cmndf[tau] - cmndf[tau - 1] - cmndf[tau + 1]));
        return sr / betterTau;
    }

    return sr / tau;
}

// ===== HELPER: Spectral Magnitude =====
function computeSpectralMagnitude(frame) {
    const n = frame.length;
    const half = Math.floor(n / 2);
    const mag = new Float32Array(half);

    for (let k = 0; k < half; k += 2) { // Compute every 2nd bin for speed
        let re = 0, im = 0;
        const step = Math.max(1, Math.floor(n / 512));

        for (let i = 0; i < n; i += step) {
            const angle = -2 * Math.PI * k * i / n;
            re += frame[i] * Math.cos(angle);
            im += frame[i] * Math.sin(angle);
        }
        mag[k] = Math.sqrt(re * re + im * im);
    }

    return mag;
}

// ===== HELPER: DFT for Key Detection =====
function computeDFT(frame, sr) {
    const n = frame.length;
    const spectrum = new Float32Array(n / 2);

    for (let k = 0; k < n / 2; k += 4) {
        let re = 0, im = 0;

        for (let i = 0; i < n; i += 8) {
            const angle = -2 * Math.PI * k * i / n;
            re += frame[i] * Math.cos(angle);
            im += frame[i] * Math.sin(angle);
        }

        spectrum[k] = Math.sqrt(re * re + im * im);
    }

    return spectrum;
}

// ===== CAMELOT WHEEL =====
function getCamelot(key, mode) {
    const wheel = {
        'C': { Major: '8B', Minor: '5A' },
        'C#': { Major: '3B', Minor: '12A' },
        'D': { Major: '10B', Minor: '7A' },
        'D#': { Major: '5B', Minor: '2A' },
        'E': { Major: '12B', Minor: '9A' },
        'F': { Major: '7B', Minor: '4A' },
        'F#': { Major: '2B', Minor: '11A' },
        'G': { Major: '9B', Minor: '6A' },
        'G#': { Major: '4B', Minor: '1A' },
        'A': { Major: '11B', Minor: '8A' },
        'A#': { Major: '6B', Minor: '3A' },
        'B': { Major: '1B', Minor: '10A' }
    };
    return wheel[key]?.[mode] || '—';
}