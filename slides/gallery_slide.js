// Week 4 gallery-walk slide: one plot, one number with its error, one sentence of surprise.
// Run:  python slides/gallery_data.py && node slides/gallery_slide.js
const path = require("path");
const fs = require("fs");
const pptxgen = require("pptxgenjs");
const { applyTheme } = require(process.env.PPTX_SKILL_DIR
  ? path.join(process.env.PPTX_SKILL_DIR, "scripts/apply_theme.js")
  : "./apply_theme.js");

const data = JSON.parse(fs.readFileSync(path.join(__dirname, "gallery_data.json"), "utf8"));
const OUT = path.join(__dirname, "gallery_slide_T3.pptx");

const THEME = {
  name: "P1-C Tissue",
  headFontFace: "Cambria",
  bodyFontFace: "Calibri",
  colors: {
    dk1: "22252B", lt1: "FFFFFF", dk2: "5A616B", lt2: "F3F1F2",
    accent1: "8E2A3B", accent2: "6B7785", accent3: "C9A227", accent4: "3F6F73",
    accent5: "B5BCC4", accent6: "D9534F", hlink: "3F6F73", folHlink: "6B7785",
  },
};
const HEX = { text: "22252B", muted: "6B7785", grid: "E3E6EA", accent: "8E2A3B", grey: "9AA3AD" };

const pres = new pptxgen();
pres.layout = "LAYOUT_16x9"; // 10 x 5.625 in
pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
pres.title = "P1-C Group 3 gallery slide";
pres.author = "Group 3";
const C = pres.SchemeColor;

pres.defineSlideMaster({
  title: "TITLE_ONLY_PH",
  background: { color: "FFFFFF" },
  objects: [
    { text: { text: "P1-C · Group 3 (T3) · Sprint 1 · analytical multilayer script",
      options: { x: 0.5, y: 5.22, w: 6.5, h: 0.28, fontSize: 10, color: C.text2, isTextBox: true, margin: 0 } } },
    { placeholder: { options: { name: "title", type: "title", x: 0.5, y: 0.3, w: 9.0, h: 0.75,
      fontSize: 26, bold: true, color: C.text1, align: "left", valign: "middle", margin: 0 }, text: "" } },
  ],
  slideNumber: { x: 9.1, y: 5.22, w: 0.4, h: 0.28, fontSize: 10, color: C.text2, align: "right" },
});

pres.addSection({ title: "Gallery walk" });
const s = pres.addSlide({ masterName: "TITLE_ONLY_PH", sectionTitle: "Gallery walk" });
s.addText("Our script matches Fresnel at all four frequencies", { placeholder: "title" });

// Plot: |Gamma| vs frequency, bare muscle (Fresnel vs script) and the skin/fat/muscle stack
const labels = data.rows.map((r) => `${r.f_MHz} MHz`);
const series = (name, key, digits, opts) => ({ type: pres.charts.LINE,
  data: [{ name, labels, values: data.rows.map((r) => +r[key].toFixed(digits)) }], options: opts });
s.addChart([
  series("Fresnel, air→muscle", "fresnel", 4, { chartColors: [HEX.grey], lineSize: 8, lineDataSymbol: "none" }),
  series("Script, air→muscle", "script", 4, { chartColors: [HEX.accent], lineSize: 0, lineDataSymbol: "circle", lineDataSymbolSize: 8 }),
  series("Script, skin 1.5 / fat 8 / muscle", "stack", 3, { chartColors: ["3F6F73"], lineSize: 2, lineDataSymbol: "circle", lineDataSymbolSize: 6 }),
], {
  x: 0.5, y: 1.2, w: 5.6, h: 3.35, objectName: "gamma-vs-frequency",
  showTitle: true, title: "|Γ| at normal incidence", titleFontSize: 12, titleColor: HEX.text, titleFontFace: "+mn-lt",
  valAxisMinVal: 0.4, valAxisMaxVal: 0.9, valAxisMajorUnit: 0.1, valAxisLabelFormatCode: "0.0",
  valAxisLabelColor: HEX.muted, catAxisLabelColor: HEX.muted,
  valAxisLabelFontSize: 10, catAxisLabelFontSize: 10, valAxisLabelFontFace: "+mn-lt", catAxisLabelFontFace: "+mn-lt",
  valGridLine: { color: HEX.grid, size: 0.75 }, catGridLine: { style: "none" },
  showLegend: true, legendPos: "b", legendFontSize: 10, legendColor: HEX.text, legendFontFace: "+mn-lt",
});

// One number with its error
s.addText("0.7625", { x: 6.5, y: 1.2, w: 3.0, h: 0.8, fontSize: 48, bold: true, color: C.accent1,
  fontFace: THEME.headFontFace, isTextBox: true, margin: 0, objectName: "stat-gamma" });
s.addText([
  { text: "|Γ| air→muscle at 2.45 GHz", options: { bold: true, breakLine: true } },
  { text: `Fresnel error ${data.gamma_err < 1e-12 ? "< 1e-12" : data.gamma_err.toExponential(1)} (slide: 0.762, 41.9 % enters)` },
], { x: 6.5, y: 2.0, w: 3.0, h: 0.6, fontSize: 12, color: C.text1, isTextBox: true, margin: 0, valign: "top" });

s.addText([
  { text: "λ/4 match test  ", options: { bold: true } },
  { text: `|Γ| = ${data.qw_gamma.toExponential(0)}`, options: { breakLine: true } },
  { text: "Power balance  ", options: { bold: true } },
  { text: `error ${data.conservation.toExponential(0)}`, options: { breakLine: true } },
  { text: "Muscle δ  ", options: { bold: true } },
  { text: `${data.delta_mm.toFixed(2)} mm (slide: 2.23 cm)` },
], { x: 6.5, y: 2.85, w: 3.0, h: 1.0, fontSize: 12, color: C.text1, isTextBox: true, margin: 0,
  valign: "top", paraSpaceAfter: 4, objectName: "other-checks" });

// One sentence of surprise
s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 6.5, y: 3.95, w: 3.0, h: 0.6, rectRadius: 0.08,
  fill: { color: C.background2 }, line: { color: C.background2 }, objectName: "surprise-card" });
s.addText([
  { text: "Surprise: ", options: { bold: true, color: C.accent1 } },
  { text: `the Week 1 δ fit on a 60 mm slab is off by ${(100 * data.fit5_err).toFixed(1)} % even on the exact field (back-face echo).` },
], { x: 6.6, y: 3.97, w: 2.8, h: 0.56, fontSize: 10, color: C.text1, isTextBox: true, margin: 0, valign: "middle" });

s.addText(`Skipping 1.5δ at the bottom brings it to ${(100 * data.fitsafe_err).toFixed(1)} %. Feko Model B rows go in the next version of this slide.`,
  { x: 0.5, y: 4.7, w: 9.0, h: 0.35, fontSize: 11, italic: true, color: C.text2, isTextBox: true, margin: 0 });

s.addNotes("Reference: Fresnel with Gabriel 1996 muscle. Test 2 is a lossless quarter-wave layer (eps 3 on eps 9). " +
  "The stack curve shows why the layered case matters: |Gamma| drops to 0.52 at 2.45 GHz.");

(async () => {
  await pres.writeFile({ fileName: OUT });
  await applyTheme(OUT, THEME);
  console.log("Wrote", OUT);
})();
