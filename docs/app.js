(() => {
  // src/chirpy-rx.js
  var ToneStencil = class {
    constructor(freq, sampleRate2, fftSize2) {
      this.freq = freq;
      this.bins = getBins(freq, sampleRate2, fftSize2, true);
    }
  };
  function getBins(freq, sampleRate2, fftSize2, multiple = false) {
    const bandwidth = sampleRate2 / fftSize2;
    let midIx = -1;
    for (let i = 0; i < fftSize2 / 2; ++i) {
      if (freq > i * bandwidth && freq <= (i + 1) * bandwidth) {
        midIx = i;
        break;
      }
    }
    if (multiple)
      return [midIx - 1, midIx, midIx + 1];
    else
      return [midIx];
  }
  var Demodulator = class {
    constructor({ sampleRate: sampleRate2, fftSize: fftSize2, toneRate: toneRate2, baseFreq: baseFreq2, freqStep: freqStep2, nFreqs: nFreqs2 }) {
      const bitSize = Math.log(nFreqs2 - 1) / Math.log(2);
      if (bitSize != Math.round(bitSize))
        throw "nFreqs must be 2^x+1, e.g., 5, 9 or 17";
      this.bitSize = bitSize;
      this.sampleRate = sampleRate2;
      this.fftSize = fftSize2;
      this.toneRate = toneRate2;
      this.sampleLenMsec = this.fftSize / this.sampleRate * 1e3;
      this.toneLenMsec = 1e3 / this.toneRate;
      this.symFreqs = [];
      for (let i = 0; i < nFreqs2; ++i)
        this.symFreqs.push(baseFreq2 + freqStep2 * i);
      this.stencils = [];
      for (const freq of this.symFreqs)
        this.stencils.push(new ToneStencil(freq, sampleRate2, fftSize2));
      console.log("Demodulator params:", { sampleRate: sampleRate2, fftSize: fftSize2, toneRate: toneRate2, baseFreq: baseFreq2, freqStep: freqStep2, nFreqs: nFreqs2 });
      const bandwidth = sampleRate2 / fftSize2;
      console.log("Bandwidth:", bandwidth);
      this.symFreqs.forEach((freq, index) => {
        let midIx = -1;
        for (let i = 0; i < fftSize2 / 2; ++i) {
          if (freq > i * bandwidth && freq <= (i + 1) * bandwidth) {
            midIx = i;
            break;
          }
        }
        console.log(`Freq ${freq} (index ${index}): midIx = ${midIx}, bins = [${midIx - 1}, ${midIx}, ${midIx + 1}]`);
      });
    }
    detecToneAt(spectra2, msec) {
      const ixAt = Math.round(msec / this.sampleLenMsec);
      const tone0 = detectTone(spectra2[ixAt - 1], this.stencils);
      const tone1 = detectTone(spectra2[ixAt], this.stencils);
      const tone2 = detectTone(spectra2[ixAt + 1], this.stencils);
      if (tone0 == tone1 || tone0 == tone2)
        return tone0;
      if (tone1 == tone2)
        return tone1;
      return -1;
    }
    findStartMsec(spectra2) {
      let firstMatchIx = -1, lastMatchIx = -1;
      for (let ix0 = 0; ix0 < spectra2.length; ++ix0) {
        const msec0 = ix0 * this.sampleLenMsec;
        const ix1 = Math.round((msec0 + this.toneLenMsec) / this.sampleLenMsec);
        const ix2 = Math.round((msec0 + 2 * this.toneLenMsec) / this.sampleLenMsec);
        const ix3 = Math.round((msec0 + 3 * this.toneLenMsec) / this.sampleLenMsec);
        if (ix3 > spectra2.length - 1)
          break;
        const tone0 = detectTone(spectra2[ix0], this.stencils);
        const tone1 = detectTone(spectra2[ix1], this.stencils);
        const tone2 = detectTone(spectra2[ix2], this.stencils);
        const tone3 = detectTone(spectra2[ix3], this.stencils);
        if (tone0 == this.symFreqs.length - 1 && tone1 == 0 && tone2 == this.symFreqs.length - 1 && tone3 == 0) {
          if (firstMatchIx == -1) {
            firstMatchIx = lastMatchIx = ix0;
          } else
            lastMatchIx = ix0;
        } else if (firstMatchIx != -1)
          break;
      }
      if (firstMatchIx == -1)
        return -1;
      const midMatchIx = Math.round((firstMatchIx + lastMatchIx) / 2);
      return Math.floor(midMatchIx * this.sampleLenMsec);
    }
  };
  var v1 = null;
  function detectTone(spectrum, stencils) {
    if (!v1 || v1.length != stencils.length)
      v1 = new Float32Array(stencils.length);
    for (let i = 0; i < v1.length; ++i)
      v1[i] = 0;
    for (let toneIx = 0; toneIx < stencils.length; ++toneIx) {
      const stencil = stencils[toneIx];
      for (const binIx of stencil.bins)
        v1[toneIx] += spectrum[binIx];
    }
    let maxVal = Number.MIN_VALUE, maxIx = -1;
    for (let i = 0; i < v1.length; ++i) {
      if (v1[i] > maxVal) {
        maxVal = v1[i];
        maxIx = i;
      }
    }
    let restSum = 0;
    for (let i = 0; i < v1.length; ++i) {
      if (i != maxIx)
        restSum += v1[i];
    }
    let ratio = maxVal / restSum;
    if (ratio >= 0.1)
      return maxIx;
    else
      return -1;
  }
  var Block = class {
    constructor(startTonePos, nTones, bytes, crc) {
      this.startTonePos = startTonePos;
      this.nTones = nTones;
      this.bytes = bytes;
      this.ascii = getAscii(bytes);
      this.crc = crc;
      this.valid = crc == getCRC8(bytes);
    }
  };
  var Decoder = class {
    constructor(tones2) {
      this.tones = tones2;
      this.blocks = decode(tones2);
      this.bytes = catBytes(this.blocks);
      this.ascii = catAscii(this.blocks);
      this.valid = true;
      for (const block of this.blocks)
        if (!block.valid)
          this.valid = false;
    }
  };
  function getAscii(bytes) {
    let res = "";
    for (const b of bytes) {
      res += String.fromCodePoint(b);
    }
    return res;
  }
  function catBytes(blocks) {
    const bytes = [];
    for (const block of blocks) {
      bytes.push(...block.bytes);
    }
    return bytes;
  }
  function catAscii(blocks) {
    let str = "";
    for (const block of blocks) {
      str += block.ascii;
    }
    return str;
  }
  function decode(tones2) {
    const blocks = [];
    if (tones2.length < 14)
      return blocks;
    if (tones2[0] != 8 || tones2[1] != 0 || tones2[2] != 8 || tones2[3] != 0)
      return blocks;
    let ix = 4;
    while (true) {
      const endIx = getBlockEndIx(tones2, ix);
      if (endIx == -1)
        break;
      const block = decodeBlock(tones2.slice(ix, endIx));
      block.startTonePos = ix;
      block.nTones = endIx - ix;
      blocks.push(block);
      ix = endIx;
    }
    return blocks;
  }
  var toneBits = [
    [0, 0, 0],
    [0, 0, 1],
    [0, 1, 0],
    [0, 1, 1],
    [1, 0, 0],
    [1, 0, 1],
    [1, 1, 0],
    [1, 1, 1]
  ];
  function getToneBits(tone) {
    return toneBits[tone % 8];
  }
  function decodeBlock(tones2, start, end) {
    const seq = tones2.slice(start, end);
    const bits = [];
    for (let i = 0; i < seq.length - 5; ++i)
      bits.push(...getToneBits(seq[i]));
    const crcBits = [
      ...getToneBits(seq[seq.length - 4]),
      ...getToneBits(seq[seq.length - 3]),
      ...getToneBits(seq[seq.length - 2])
    ];
    const bytes = getBytes(bits);
    const crcBytes = getBytes(crcBits);
    return new Block(start, end - start, bytes, crcBytes[0]);
  }
  function getBytes(bits) {
    const res = [];
    for (let i = 0; i + 8 <= bits.length; i += 8) {
      let val = 0;
      for (let j = 0; j < 8; ++j) {
        val <<= 1;
        val += bits[i + j];
      }
      res.push(val);
    }
    return res;
  }
  function getBlockEndIx(tones2, startIx) {
    for (let i = startIx + 4; i < tones2.length; ++i) {
      if (tones2[i] == 8 && tones2[i - 4] == 8) {
        return i + 1;
      }
    }
    return -1;
  }
  function getCRC8(bytes) {
    let crc = 0;
    for (const b of bytes)
      crc = updateCRC(b, crc);
    return crc;
    function updateCRC(nextByte, crc2) {
      for (let j = 0; j < 8; j++) {
        let mix = (crc2 ^ nextByte) & 1;
        crc2 >>= 1;
        if (mix)
          crc2 ^= 140;
        nextByte >>= 1;
      }
      return crc2;
    }
  }

  // src/wav-encoder.js
  var min = Math.min;
  var max = Math.max;
  function setString(view, offset, str) {
    let len = str.length;
    for (let i = 0; i < len; ++i)
      view.setUint8(offset + i, str.charCodeAt(i));
  }
  var WAVEncoder = class {
    constructor(sampleRate2, numChannels) {
      this.sampleRate = sampleRate2;
      this.numChannels = numChannels;
      this.numSamples = 0;
      this.dataViews = [];
    }
    encode(buffer) {
      let len = buffer[0].length, nCh = this.numChannels, view = new DataView(new ArrayBuffer(len * nCh * 2)), offset = 0;
      for (let i = 0; i < len; ++i)
        for (let ch = 0; ch < nCh; ++ch) {
          let x = buffer[ch][i] * 32767;
          view.setInt16(offset, x < 0 ? max(x, -32768) : min(x, 32767), true);
          offset += 2;
        }
      this.dataViews.push(view);
      this.numSamples += len;
    }
    finish(mimeType) {
      let dataSize = this.numChannels * this.numSamples * 2, view = new DataView(new ArrayBuffer(44));
      setString(view, 0, "RIFF");
      view.setUint32(4, 36 + dataSize, true);
      setString(view, 8, "WAVE");
      setString(view, 12, "fmt ");
      view.setUint32(16, 16, true);
      view.setUint16(20, 1, true);
      view.setUint16(22, this.numChannels, true);
      view.setUint32(24, this.sampleRate, true);
      view.setUint32(28, this.sampleRate * 4, true);
      view.setUint16(32, this.numChannels * 2, true);
      view.setUint16(34, 16, true);
      setString(view, 36, "data");
      view.setUint32(40, dataSize, true);
      this.dataViews.unshift(view);
      let blob = new Blob(this.dataViews, { type: "audio/wav" });
      this.cleanup();
      return blob;
    }
    cancel() {
      this.cleanup();
    }
    cleanup() {
      delete this.dataViews;
    }
    canEncode() {
      return this.hasOwnProperty("dataViews");
    }
  };

  // src/fft.js
  function FourierTransform(bufferSize, sampleRate2) {
    this.bufferSize = bufferSize;
    this.sampleRate = sampleRate2;
    this.bandwidth = 2 / bufferSize * sampleRate2 / 2;
    this.spectrum = new Float64Array(bufferSize / 2);
    this.real = new Float64Array(bufferSize);
    this.imag = new Float64Array(bufferSize);
    this.peakBand = 0;
    this.peak = 0;
    this.getBandFrequency = function(index) {
      return this.bandwidth * index + this.bandwidth / 2;
    };
    this.calculateSpectrum = function() {
      let spectrum = this.spectrum, real = this.real, imag = this.imag, bSi = 2 / this.bufferSize, sqrt = Math.sqrt, rval, ival, mag;
      for (let i = 0, N = bufferSize / 2; i < N; i++) {
        rval = real[i];
        ival = imag[i];
        mag = bSi * sqrt(rval * rval + ival * ival);
        if (mag > this.peak) {
          this.peakBand = i;
          this.peak = mag;
        }
        spectrum[i] = mag;
      }
    };
  }
  var FFT = class {
    constructor(bufferSize, sampleRate2) {
      FourierTransform.call(this, bufferSize, sampleRate2);
      this.reverseTable = new Uint32Array(bufferSize);
      let limit = 1;
      let bit = bufferSize >> 1;
      let i;
      while (limit < bufferSize) {
        for (i = 0; i < limit; i++) {
          this.reverseTable[i + limit] = this.reverseTable[i] + bit;
        }
        limit = limit << 1;
        bit = bit >> 1;
      }
      this.sinTable = new Float64Array(bufferSize);
      this.cosTable = new Float64Array(bufferSize);
      for (i = 0; i < bufferSize; i++) {
        this.sinTable[i] = Math.sin(-Math.PI / i);
        this.cosTable[i] = Math.cos(-Math.PI / i);
      }
    }
    forward(buffer) {
      let bufferSize = this.bufferSize, cosTable = this.cosTable, sinTable = this.sinTable, reverseTable = this.reverseTable, real = this.real, imag = this.imag, spectrum = this.spectrum;
      let k = Math.floor(Math.log(bufferSize) / Math.LN2);
      if (Math.pow(2, k) !== bufferSize) {
        throw "Invalid buffer size, must be a power of 2.";
      }
      if (bufferSize !== buffer.length) {
        throw "Supplied buffer is not the same size as defined FFT. FFT Size: " + bufferSize + " Buffer Size: " + buffer.length;
      }
      let halfSize = 1, phaseShiftStepReal, phaseShiftStepImag, currentPhaseShiftReal, currentPhaseShiftImag, off, tr, ti, tmpReal, i;
      for (i = 0; i < bufferSize; i++) {
        real[i] = buffer[reverseTable[i]];
        imag[i] = 0;
      }
      while (halfSize < bufferSize) {
        phaseShiftStepReal = cosTable[halfSize];
        phaseShiftStepImag = sinTable[halfSize];
        currentPhaseShiftReal = 1;
        currentPhaseShiftImag = 0;
        for (let fftStep = 0; fftStep < halfSize; fftStep++) {
          i = fftStep;
          while (i < bufferSize) {
            off = i + halfSize;
            tr = currentPhaseShiftReal * real[off] - currentPhaseShiftImag * imag[off];
            ti = currentPhaseShiftReal * imag[off] + currentPhaseShiftImag * real[off];
            real[off] = real[i] - tr;
            imag[off] = imag[i] - ti;
            real[i] += tr;
            imag[i] += ti;
            i += halfSize << 1;
          }
          tmpReal = currentPhaseShiftReal;
          currentPhaseShiftReal = tmpReal * phaseShiftStepReal - currentPhaseShiftImag * phaseShiftStepImag;
          currentPhaseShiftImag = tmpReal * phaseShiftStepImag + currentPhaseShiftImag * phaseShiftStepReal;
        }
        halfSize = halfSize << 1;
      }
      return this.calculateSpectrum();
    }
    inverse(real, imag) {
      let bufferSize = this.bufferSize, cosTable = this.cosTable, sinTable = this.sinTable, reverseTable = this.reverseTable, spectrum = this.spectrum;
      real = real || this.real;
      imag = imag || this.imag;
      let halfSize = 1, phaseShiftStepReal, phaseShiftStepImag, currentPhaseShiftReal, currentPhaseShiftImag, off, tr, ti, tmpReal, i;
      for (i = 0; i < bufferSize; i++) {
        imag[i] *= -1;
      }
      let revReal = new Float64Array(bufferSize);
      let revImag = new Float64Array(bufferSize);
      for (i = 0; i < real.length; i++) {
        revReal[i] = real[reverseTable[i]];
        revImag[i] = imag[reverseTable[i]];
      }
      real = revReal;
      imag = revImag;
      while (halfSize < bufferSize) {
        phaseShiftStepReal = cosTable[halfSize];
        phaseShiftStepImag = sinTable[halfSize];
        currentPhaseShiftReal = 1;
        currentPhaseShiftImag = 0;
        for (let fftStep = 0; fftStep < halfSize; fftStep++) {
          i = fftStep;
          while (i < bufferSize) {
            off = i + halfSize;
            tr = currentPhaseShiftReal * real[off] - currentPhaseShiftImag * imag[off];
            ti = currentPhaseShiftReal * imag[off] + currentPhaseShiftImag * real[off];
            real[off] = real[i] - tr;
            imag[off] = imag[i] - ti;
            real[i] += tr;
            imag[i] += ti;
            i += halfSize << 1;
          }
          tmpReal = currentPhaseShiftReal;
          currentPhaseShiftReal = tmpReal * phaseShiftStepReal - currentPhaseShiftImag * phaseShiftStepImag;
          currentPhaseShiftImag = tmpReal * phaseShiftStepImag + currentPhaseShiftImag * phaseShiftStepReal;
        }
        halfSize = halfSize << 1;
      }
      let buffer = new Float64Array(bufferSize);
      for (i = 0; i < bufferSize; i++) {
        buffer[i] = real[i] / bufferSize;
      }
      return buffer;
    }
  };

  // src/base64.js
  function toBase64(bytes) {
    let base64 = "";
    let encodings = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let byteLength = bytes.length;
    let byteRemainder = byteLength % 3;
    let mainLength = byteLength - byteRemainder;
    let a, b, c, d;
    let chunk;
    for (let i = 0; i < mainLength; i = i + 3) {
      chunk = bytes[i] << 16 | bytes[i + 1] << 8 | bytes[i + 2];
      a = (chunk & 16515072) >> 18;
      b = (chunk & 258048) >> 12;
      c = (chunk & 4032) >> 6;
      d = chunk & 63;
      base64 += encodings[a] + encodings[b] + encodings[c] + encodings[d];
    }
    if (byteRemainder == 1) {
      chunk = bytes[mainLength];
      a = (chunk & 252) >> 2;
      b = (chunk & 3) << 4;
      base64 += encodings[a] + encodings[b] + "==";
    } else if (byteRemainder == 2) {
      chunk = bytes[mainLength] << 8 | bytes[mainLength + 1];
      a = (chunk & 64512) >> 10;
      b = (chunk & 1008) >> 4;
      c = (chunk & 15) << 2;
      base64 += encodings[a] + encodings[b] + encodings[c] + "=";
    }
    let res = "";
    for (let i = 0; i < base64.length; ++i) {
      if (i > 0 && i % 76 == 0)
        res += "\n";
      res += base64[i];
    }
    return res;
  }

  // src/content-shared.js
  function decodeDateTime(bytes) {
    let year = (bytes[0] & 252) >> 2;
    year += 2020;
    let month = ((bytes[0] & 3) << 2) + ((bytes[1] & 192) >> 6);
    let day = (bytes[1] & 62) >> 1;
    let hour = ((bytes[1] & 1) << 4) + ((bytes[2] & 240) >> 4);
    let minute = ((bytes[2] & 15) << 2) + ((bytes[3] & 192) >> 6);
    let second = bytes[3] & 63;
    return new Date(year, month - 1, day, hour, minute, second);
  }
  function dateToIso(date) {
    let res = date.getFullYear() + "-" + (date.getMonth() + 1).toString().padStart(2, "0") + "-";
    res += date.getDate().toString().padStart(2, "0") + "T";
    res += date.getHours().toString().padStart(2, "0") + ":";
    res += date.getMinutes().toString().padStart(2, "0") + ":";
    res += date.getSeconds().toString().padStart(2, "0");
    return res;
  }

  // src/content-activity.js
  function interpretActivity(bytes, elmRes, elmLnk) {
    if (bytes.length <= 2 || bytes[0] != 39 || bytes[1] != 0)
      return false;
    const itmLen = 9;
    const items = [];
    for (let ix = 2; ix < bytes.length; ix += itmLen) {
      const itmBytes = bytes.slice(ix, ix + itmLen);
      items.push(new DecodedActivity(itmBytes));
    }
    let tableTxt = "time	activity_code	activity_name	total_duration	pause_duration\n";
    for (const itm of items) {
      tableTxt += dateToIso(itm.start) + "	" + itm.typeNum + "	" + itm.type + "	";
      tableTxt += secsToDuration(itm.totalSec) + "	" + secsToDuration(itm.pauseSec);
      tableTxt += "\n";
    }
    elmRes.innerHTML = "<pre></pre>";
    elmRes.querySelector("pre").innerText = tableTxt;
    elmLnk.innerText = "Activity";
    return true;
  }
  var DecodedActivity = class {
    constructor(bytes) {
      this.start = decodeDateTime(bytes);
      this.totalSec = bytes[4] * 256 + bytes[5];
      this.pauseSec = bytes[6] * 256 + bytes[7];
      this.typeNum = bytes[8];
      this.type = "UNKNOWN";
      if (this.typeNum < activities.length)
        this.type = activities[this.typeNum];
    }
  };
  var activities = [
    "BIKE",
    "WALK",
    "RUN",
    "DANCE",
    "YOGA",
    "CROSSFIT",
    "SWIM",
    "ELLIPTICAL",
    "GYM",
    "ROW",
    "SOCCER",
    "FOOTBALL",
    "BALLGAME",
    "SKI"
  ];
  function secsToDuration(val) {
    let secs = val % 60;
    val -= secs;
    val /= 60;
    let mins = val % 60;
    val -= mins;
    let hours = val / 60;
    let res = mins.toString().padStart(2, "0") + ":" + secs.toString().padStart(2, "0");
    if (hours == 0)
      return res;
    res = hours.toString() + ":" + res;
    return res;
  }

  // src/content-nanosec.js
  function interpretNanosecIni(bytes, elmRes, elmLnk) {
    if (bytes.length < 18 + 2 || bytes[0] != 192 || bytes[1] != 0)
      return false;
    const info = new DecodedNanosecIni(bytes.slice(2, bytes.length));
    let txt = "";
    txt += "Correction profile:     " + info.correctionProfile + " (" + info.correctionPorfileStr + ")\n";
    txt += "Frequency correction:   " + (info.freqCorrection / 100).toFixed(2) + "\n";
    txt += "Center temperature:     " + (info.centerTemperature / 100).toFixed(2) + "\n";
    txt += "Quadratic temp coeff:  -" + info.quadraticTempCo + "e-5\n";
    txt += "Cubic temp coeff:       " + info.cubicTempCo + "e-7\n";
    txt += "Correction cadence:     " + info.correctionCadence + "\n";
    txt += "Last correction time:   " + dateToIso(info.lastCorrectionTime) + "\n";
    txt += "Annual aging:           " + (info.annualAgingPPA / 100).toFixed(2) + " PPA\n";
    elmRes.innerHTML = "<pre></pre>";
    elmRes.querySelector("pre").innerText = txt;
    elmLnk.innerText = "Nanosec.ini";
    return true;
  }
  var DecodedNanosecIni = class {
    constructor(bytes) {
      this.correctionProfile = bytes[0];
      this.correctionPorfileStr = profiles[this.correctionProfile];
      this.freqCorrection = (bytes[3] << 8) + bytes[2];
      this.centerTemperature = (bytes[5] << 8) + bytes[4];
      this.quadraticTempCo = (bytes[7] << 8) + bytes[6];
      this.cubicTempCo = (bytes[9] << 8) + bytes[8];
      this.correctionCadence = bytes[10];
      const lastCurrTimeUnix = (bytes[15] << 24) + (bytes[14] << 16) + (bytes[13] << 8) + bytes[12];
      this.lastCorrectionTime = new Date(lastCurrTimeUnix * 1e3);
      this.annualAgingPPA = (bytes[17] << 8) + bytes[16];
    }
  };
  var profiles = [
    "static hardware correction",
    "static correction with dithering",
    "datasheet quadratic correction",
    "cubic correction conservative",
    "cubic correction finetuned"
  ];

  // src/content-ascii.js
  function interpretAscii(bytes, elmRes, elmLnk) {
    let canBeAscii = true;
    for (const b of bytes)
      if (b == 0 || b > 127)
        canBeAscii = false;
    if (!canBeAscii)
      return false;
    let ascii = "";
    for (const b of bytes) {
      ascii += String.fromCodePoint(b);
    }
    elmRes.innerHTML = "<pre></pre>";
    elmRes.querySelector("pre").innerText = ascii;
    elmLnk.innerText = "ASCII";
    return true;
  }

  // src/content.js
  var parserFuns = [
    interpretActivity,
    interpretNanosecIni,
    interpretAscii
  ];
  function interpretContent(bytes, elmRes, elmLnk) {
    for (const fun of parserFuns) {
      if (fun(bytes, elmRes, elmLnk))
        return true;
    }
    return false;
  }

  // src/app.js
  var showTest = false;
  var testFileName = "data-06.wav";
  var gainVal = 10;
  var toneRate = 64 / 3;
  var baseFreq = 2500;
  var freqStep = 250;
  var nFreqs = 9;
  var fftSize = 512;
  var audioCtx;
  var stream;
  var gain;
  var source;
  var wavEncoder;
  var mediaRecorder;
  var recordedChunks;
  var chunks;
  var nSamples;
  var sampleRate;
  var spectra;
  var demodulator;
  var startMsec;
  var tones;
  var decoder;
  var elms = {
    btnAudio: null,
    ctrlAudioTop: null,
    lblLength: null,
    lnkWav: null,
    lnkTest: null,
    ctrlDecoding: null,
    decodingStatus: null,
    decodingRes: null,
    lnkTones: null,
    lnkBlocks: null,
    lnkBase64: null,
    lnkContent: null,
    resTones: null,
    resBlocks: null,
    resBase64: null,
    resContent: null
  };
  document.addEventListener("DOMContentLoaded", () => {
    initUI();
  });
  function initUI() {
    for (const key in elms)
      elms[key] = document.getElementById(key);
    elms.btnAudio.addEventListener("click", onBtnAudioClick);
    if (showTest) {
      elms.lnkTest.classList.add("visible");
      elms.lnkTest.addEventListener("click", onTestClick);
    }
    elms.lnkTones.addEventListener("click", () => setCtrlDecodingTab("tones"));
    elms.lnkBlocks.addEventListener("click", () => setCtrlDecodingTab("blocks"));
    elms.lnkBase64.addEventListener("click", () => setCtrlDecodingTab("base64"));
    elms.lnkContent.addEventListener("click", () => setCtrlDecodingTab("content"));
  }
  function setCtrlDecodingTab(tab) {
    ["tones", "blocks", "base64", "content"].forEach((cls) => elms.ctrlDecoding.classList.remove(cls));
    elms.ctrlDecoding.classList.add(tab);
  }
  function onTestClick() {
    audioCtx = new (window.AudioContext || window.webkitAudioContext)();
    gain = audioCtx.createGain();
    gain.connect(audioCtx.destination);
    gain.gain.setValueAtTime(gainVal, audioCtx.currentTime);
    const req = new XMLHttpRequest();
    req.open("GET", testFileName, true);
    req.responseType = "arraybuffer";
    req.onload = function() {
      const audioData = req.response;
      audioCtx.decodeAudioData(audioData).then((buf) => {
        const data = buf.getChannelData(0);
        chunks = [];
        sampleRate = buf.sampleRate;
        nSamples = data.length;
        let pos = 0;
        while (pos < data.length) {
          const chunkSize = Math.min(4096, data.length - pos);
          const chunk = new Float32Array(chunkSize);
          for (let i = 0; i < chunkSize; ++i) {
            chunk[i] = data[pos + i];
          }
          chunks.push(chunk);
          pos += chunkSize;
        }
        setAutdioBtnClass("done");
        elms.lnkTest.classList.remove("visible");
        elms.ctrlAudioTop.innerText = "Loaded audio file";
        startProcessing();
      }).catch((err) => {
        alert("Error decoding audio data: " + err.err);
      });
    };
    req.send();
  }
  function setAutdioBtnClass(cls) {
    ["enable", "record", "stop", "done"].forEach((cls2) => elms.btnAudio.classList.remove(cls2));
    elms.btnAudio.classList.add(cls);
  }
  function onBtnAudioClick() {
    if (elms.btnAudio.classList.contains("enable")) {
      navigator.mediaDevices.getUserMedia({ audio: true }).then((s) => {
        audioCtx = new (window.AudioContext || window.webkitAudioContext)();
        stream = s;
        setAutdioBtnClass("record");
        elms.ctrlAudioTop.innerText = "Press to record transmission";
      }).catch((err) => {
        elms.ctrlAudioTop.innerText = "Failed to enable microphone ;-(";
      });
    } else if (elms.btnAudio.classList.contains("record")) {
      let updateTime = function() {
        const elapsed = new Date() - startTime;
        let secs = Math.floor(elapsed / 1e3);
        let mins = Math.floor(secs / 60);
        secs -= mins * 60;
        elms.lblLength.innerText = mins + ":" + secs.toString().padStart(2, "0");
        if (elms.btnAudio.classList.contains("stop")) {
          setTimeout(updateTime, 50);
        }
      };
      elms.lnkTest.classList.remove("visible");
      source = audioCtx.createMediaStreamSource(stream);
      gain = audioCtx.createGain();
      source.connect(gain);
      recordedChunks = [];
      mediaRecorder = new MediaRecorder(stream);
      mediaRecorder.ondataavailable = function(event) {
        if (event.data.size > 0) {
          recordedChunks.push(event.data);
        }
      };
      mediaRecorder.onstop = function() {
        const webmBlob = new Blob(recordedChunks, { "type": "audio/webm" });
        webmBlob.arrayBuffer().then((arrayBuffer) => {
          audioCtx.decodeAudioData(arrayBuffer).then((audioBuffer) => {
            const data = audioBuffer.getChannelData(0);
            sampleRate = audioBuffer.sampleRate;
            wavEncoder = new WAVEncoder(audioBuffer.sampleRate, 1);
            wavEncoder.encode([data]);
            const wavBlob = wavEncoder.finish();
            const url = window.URL.createObjectURL(wavBlob);
            elms.lnkWav.href = url;
            elms.lnkWav.download = "chirpy.wav";
            elms.lnkWav.style.display = "inline";
            chunks = [];
            let pos = 0;
            while (pos < data.length) {
              const chunkSize = Math.min(4096, data.length - pos);
              const chunk = new Float32Array(chunkSize);
              for (let i = 0; i < chunkSize; ++i) {
                chunk[i] = data[pos + i];
              }
              chunks.push(chunk);
              pos += chunkSize;
            }
            startProcessing();
          }).catch((err) => {
            console.error("Error decoding recorded audio data:", err);
            elms.decodingStatus.innerText = "Failed to process recording.";
          });
        });
      };
      mediaRecorder.start();
      setAutdioBtnClass("stop");
      elms.ctrlAudioTop.innerText = "Press to finish recording";
      const startTime = new Date();
      updateTime();
    } else if (elms.btnAudio.classList.contains("stop")) {
      setAutdioBtnClass("done");
      elms.ctrlAudioTop.innerText = "Recording finished";
      source.disconnect();
      mediaRecorder.stop();
    }
  }
  function startProcessing() {
    elms.ctrlDecoding.classList.add("visible");
    const fft = new FFT(fftSize, sampleRate);
    spectra = [];
    let chunkIx = 0, posInChunk = 0;
    const frame = new Float32Array(fftSize);
    const framesPerIter = 1e3;
    decodeSome();
    function decodeSome() {
      let dataOver = false;
      for (let fc = 0; fc < framesPerIter && !dataOver; ++fc) {
        for (let i = 0; i < fftSize; ++i) {
          if (posInChunk == chunks[chunkIx].length) {
            posInChunk = 0;
            chunkIx += 1;
          }
          if (chunkIx == chunks.length) {
            dataOver = true;
            break;
          }
          frame[i] = chunks[chunkIx][posInChunk] * gainVal;
          ++posInChunk;
        }
        if (dataOver)
          break;
        fft.forward(frame);
        let s = new Float32Array(fft.spectrum.length);
        for (let i = 0; i < s.length; ++i)
          s[i] = fft.spectrum[i];
        spectra.push(s);
      }
      if (dataOver)
        startDemodulating();
      else {
        setTimeout(decodeSome, 10);
      }
    }
  }
  function startDemodulating() {
    demodulator = new Demodulator({
      sampleRate,
      fftSize,
      toneRate,
      baseFreq,
      freqStep,
      nFreqs
    });
    startMsec = demodulator.findStartMsec(spectra);
    if (startMsec == -1) {
      elms.decodingStatus.innerText = "No message found.";
      elms.resTones.innerHTML += "<p>No Start-Of-Message sequence detected. Cannot decode transmission.</p>";
      return;
    }
    tones = [];
    let tonePos = 0;
    const tonesPerIter = 500;
    const recLenMsec = Math.round(nSamples / sampleRate * 1e3);
    demodulateSome();
    function demodulateSome() {
      for (let tc = 0; tc < tonesPerIter; ++tc) {
        const msec = startMsec + tonePos * demodulator.toneLenMsec;
        if (msec + 200 > recLenMsec) {
          decodeTones(startMsec, null);
          return;
        }
        const tone = demodulator.detecToneAt(spectra, msec);
        tones.push(tone);
        if (doesEndInEOM(tones, demodulator.symFreqs.length - 1)) {
          decodeTones(startMsec, msec);
          return;
        }
        ++tonePos;
      }
      setTimeout(demodulateSome, 10);
    }
    function doesEndInEOM(tones2, signalToneIx) {
      if (tones2.length < 3)
        return false;
      for (let i = 0; i < 3; ++i) {
        if (tones2[tones2.length - i - 1] != signalToneIx)
          return false;
      }
      return true;
    }
  }
  function decodeTones(startMsec2, endMsec) {
    console.log("Detected tones:", tones);
    const startSecStr = (startMsec2 / 1e3).toFixed(2);
    if (!endMsec) {
      elms.resTones.innerHTML += `<p>Start of message: ${startSecStr}<br/>No End-Of-Message sequence detected</p>`;
    } else {
      const endSecStr = (endMsec / 1e3).toFixed(2);
      elms.resTones.innerHTML += `<p>Start of message: ${startSecStr}<br/>End of message: ${endSecStr}</p>`;
    }
    let tonesStr = "";
    for (const t of tones) {
      if (tonesStr != "")
        tonesStr += " ";
      tonesStr += t;
    }
    elms.resTones.innerHTML += `<pre>${tonesStr}</pre>`;
    decoder = new Decoder(tones);
    let blocksHtml = "";
    for (let i = 0; i < decoder.blocks.length; ++i) {
      const block = decoder.blocks[i];
      if (block.valid)
        blocksHtml += `<span class='valid'>Block ${i} VALID</span>`;
      else
        blocksHtml += `<span class='invalid'>Block ${i} INVALID</span>`;
      blocksHtml += "\nTones:";
      for (let j = block.startTonePos; j < block.startTonePos + block.nTones; ++j)
        blocksHtml += " " + tones[j];
      blocksHtml += "\nBytes:";
      for (const b of block.bytes)
        blocksHtml += " 0x" + b.toString(16).padStart(2, "0");
      blocksHtml += "\nCRC: 0x" + block.crc.toString(16).padStart(2, "0") + "\n\n";
    }
    elms.resBlocks.innerHTML = `<pre>${blocksHtml}</pre>`;
    elms.lnkBlocks.classList.add("visible");
    setCtrlDecodingTab("blocks");
    if (!decoder.valid) {
      elms.decodingStatus.innerText = "Message cannot be reconstructed: invalid CRC in one or more blocks.";
      return;
    }
    const base64 = toBase64(decoder.bytes);
    elms.resBase64.innerHTML = `
<p>
  This is the transmission's data content encoded in Base64. If you're not sure how to decode it,
  you can use an online tool like
  <a href="https://base64.guru/converter/decode" target="_blank" rel="noreferrer">Base64.guru</a>.
</p>
<pre>${base64}</pre>`;
    elms.lnkBase64.classList.add("visible");
    setCtrlDecodingTab("base64");
    elms.decodingStatus.innerText = "Message successfully decoded.";
    if (interpretContent(decoder.bytes, elms.resContent, elms.lnkContent)) {
      elms.lnkContent.classList.add("visible");
      setCtrlDecodingTab("content");
    }
  }
})();
//# sourceMappingURL=app.js.map
