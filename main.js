const _id = id => document.getElementById(id);
function $(p,...args) {
  if (p===null) {
    if (args[0].constructor !== String) throw new Error('expected tag name');
    p = document.createElement(args.shift());
  }
  for (let x of args) {
    if (x.constructor === String) {
      p = p.appendChild( (p instanceof SVGElement || x==='svg')
        ? document.createElementNS('http://www.w3.org/2000/svg',x)
        : document.createElement(x)
      );
    } else if (x.constructor === Array) {
      x = x.filter(x=>!!x);
      if (x.length) p.classList.add(...x);
    } else if (x.constructor === Object) {
      for (const [key,val] of Object.entries(x)) {
        if (key==='style') {
          for (const [skey,sval] of Object.entries(val)) {
            if (sval!==null)
              p.style[skey] = sval;
            else
              p.style.removeProperty(skey);
          }
        } else {
          if (p instanceof SVGElement)
            p.setAttributeNS(null,key,val);
          else
            p.setAttribute(key,val);
        }
      }
    }
  }
  return p;
}
function clear(x) {
  for (let c; c = x.firstChild; ) x.removeChild(c);
  return x;
}
const round = x => x.toFixed(4).replace(/\.?0*$/,'');
const last = xs => xs[xs.length-1];

function* enumerate(xs,i=0) {
  for (const x of xs) yield [i++, x];
}

var lumi, nom_lumi = 140.429, var_name = 'pT_yy', table;

var em;

function resize_table() {
  const w  = table.parentNode.clientWidth;
  const wt = table.offsetWidth;
  let x = parseFloat(
    window.getComputedStyle(table,null).getPropertyValue('font-size')
  );
  if (w < wt || x < em) {
    x *= w/wt;
    if (x > em) x = em;
    table.style.fontSize = `${x}px`;
  }
}

document.addEventListener('DOMContentLoaded', () => {
  em = parseFloat(
    window.getComputedStyle(document.body,null).getPropertyValue('font-size')
  );
  table = _id('t');
  const form = _id('q');
  const H = form.querySelector('[name="H"]');
  const V = form.querySelector('[name="V"]');
  for (const x of Object.keys(vars).sort().reverse())
    $(H,'option',{'value':x}).textContent = x;
  H.onchange = function(e) {
    clear(V);
    for (const x of vars[this.value]) {
      const opt = $(V,'option',{'value':x});
      if (x===var_name) opt.selected = 'selected';
      opt.textContent = x;
    }
  };
  V.onchange = function(e) { var_name = this.value; };
  H.onchange();
  V.onchange();

  function submit(e) {
    if (e) e.preventDefault();
    let q = '', m, v;

    m = form.querySelector('[name="L"]');
    v = parseFloat(m.value);
    if (!(v > 0)) v = nom_lumi;
    lumi = v;
    m.value = v = `${v}`;
    q += '?L='+v;

    q += '&H='+H.value;
    q += '&V='+var_name;

    m = form.querySelector('[name="edges"]');
    v = m.value.split(/[\s,;:]+/).map(
      x => parseFloat(
        x.replace(/(?:\binf(?:(?:ini)?ty)?\b|∞)/gi,'Infinity')
         .replace(/−/g,'-')
      )
    ).sort((a,b)=>a-b).filter((x,i,xs) => !isNaN(x) && (i==0 || x!==xs[i-1]));
    if (v.length == 0) {
      try {
        v = default_binning[var_name][0][1];
        if (!(v instanceof Array) || v.length == 0)
          throw 'Incorrectly formatted default binning';
      } catch (err) {
        console.error(err);
        v = [ -Infinity ];
      }
    }
    if (v.length == 1) {
      v = v[0]!==Infinity ? [ v[0], Infinity ] : [ -Infinity, Infinity ];
    }
    m.value = v.join(' ').replace(/\bInfinity\b/g,'∞');
    q += '&edges=' + v.join('+').replace(/\bInfinity\b/g,'inf');

    m = form.querySelector('[name="Bf"]');
    q += '&Bf='+m.value;

    m = form.querySelector('[name="Bdeg"]');
    v = parseInt(m.value);
    if (isNaN(v)) v = 2;
    q += `&Bdeg=${v}`;

    m = form.querySelector('[name="Bdiv"]');
    v = parseInt(m.value);
    if (isNaN(v)) v = 1;
    q += `&Bdiv=${v}`;

    m = form.querySelector('[name="Sdiv"]');
    v = parseInt(m.value);
    if (isNaN(v)) v = 1;
    q += `&Sdiv=${v}`;

    m = form.querySelector('[name="nV"]');
    v = parseInt(m.value);
    if (isNaN(v)) v = 100;
    q += `&nV=${v}`;

    let u = q;

    m = form.querySelector('[name="method"]');
    const method = m.value;
    u += '&method='+method;

    m = form.querySelector('[name="click"]');
    if (m.checked) u += '&click';

    m = form.querySelector('[name="unc"]');
    if (m.checked) u += '&unc';

    console.log(q);
    console.log(u);

    clear(table);
    const cls = [
      null,null,'unc','unc',
      null,null,null,'unc','unc',
      null,null,null,null
    ];
    for (const row of [[
      '','[105,160]','[121,129]','MC unc.','cnt. unc.',
      '[105,121)','(129,160]','[121,129]','fit unc.','cnt. unc.',
      'signif','signif','','reco'
    ],[
      'bin','sig','sig','\u221a\u2211w\u00B2','\u221asig',
      'L bkg','R bkg','bkg','','\u221abkg',
      's/\u221a(s+b)','Cowan','s/(s+b)','purity'
    ]]) {
      const tr = $(table,'tr');
      for (const [i,col] of enumerate(row)) {
        const td = $(tr,'td',[cls[i]]);
        td.textContent = col;
      }
    }
    resize_table();

    fetch('req.php'+q, { method: 'GET' })
    .then(r => r.json())
    .then(r => {
      console.log(r);
      _id('time').textContent = `(${r['time'].toFixed(2)} sec)`;
      const edges = r.edges.map(x => x==='inf' ? '∞' : x==='-inf' ? '-∞' : x);

      for (const [b,bin] of enumerate(r.bins)) {
        const tr = $(table,'tr');
        $(tr,'td').textContent = `[${edges[b]},${edges[b+1]})`;
        { let sum_all=0, sum_sig=0;
          const hist = bin.S.hist,
                n3 = hist.length,
                n1 = n3*(121-105)/(160-105),
                n2 = n3*(129-105)/(160-105);
          for (const x of hist)
            sum_all += x;
          for (let i=n1; i<n2; ++i)
            sum_sig += hist[i];
          sum_all *= lumi;
          sum_sig *= lumi;
          $(tr,'td').textContent = sum_all.toFixed(2);
          $(tr,'td').textContent = sum_sig.toFixed(2);
          $(tr,'td');
          $(tr,'td').textContent = (100*Math.sqrt(sum_sig)/sum_sig).toFixed(2)+'%';
        }
        { let Lbkg=0, bkg=0, Rbkg=0;
          const hist = bin.B.hist,
                n3 = hist.length*(160-105)/(121-105+160-129),
                n1 = n3*(121-105)/(160-105),
                n2 = n3 - n3*(129-121)/(160-105);
          let i=0;
          for (; i<n1; ++i) Lbkg += hist[i];
          for (; i<n2; ++i) Rbkg += hist[i];
          $(tr,'td').textContent = Lbkg.toFixed(0);
          $(tr,'td').textContent = Rbkg.toFixed(0);
          $(tr,'td').textContent =  bkg.toFixed(0);
          $(tr,'td');
          $(tr,'td').textContent = '%';
        }
        $(tr,'td');
        $(tr,'td');
        $(tr,'td');
        $(tr,'td');
      }
      resize_table();
    });
  }
  form.addEventListener('submit',submit);
  submit();
});

window.addEventListener('resize',resize_table);
