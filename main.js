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
const last = xs => xs[xs.length-1];

function* enumerate(xs,i=0) {
  for (const x of xs) yield [i++, x];
}

var lumi, var_name = 'pT_yy', table, em;

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
function fix_edges(edges) {
  return edges.map(
    x => typeof x === 'string' ? parseFloat(
      x.replace(/(?:\binf(?:(?:ini)?ty)?\b|∞)/gi,'Infinity')
       .replace(/−/g,'-')
    ) : x
  ).sort((a,b)=>a-b).filter((x,i,xs) => !isNaN(x) && (i==0 || x!==xs[i-1]));
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

    q += '?H='+H.value;
    q += '&V='+var_name;

    m = form.querySelector('[name="edges"]');
    v = fix_edges(m.value.split(/[\s,;:]+/));
    if (v.length == 0) {
      try {
        v = default_binning[var_name][0][1];
        if (!(v instanceof Array) || v.length == 0)
          throw 'Incorrectly formatted default binning';
        v = fix_edges(v);
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

    m = form.querySelector('[name="L"]');
    lumi = parseFloat(m.value);
    if (!(lumi > 0)) {
      lumi = 0;
      m.value = '';
    } else {
      v = `${lumi}`;
      m.value = v;
      u += '&L='+v;
    }

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
      '[105,121]','[129,160]','[121,129]','fit unc.','cnt. unc.',
      'signif','signif','','reco'
    ],[
      'bin','sig','sig','\u221a\u2211w\u00B2','\u221asig',
      'L data','R data','bkg','','\u221abkg',
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
    .then(r => {
      try {
        return r.json();
      } catch (e) {
        console.log(r.text());
        throw e;
      }
    })
    .then(r => {
      console.log(r);
      _id('time').textContent = `(${r['time'].toFixed(2)} sec)`;

      const edges = r.edges.map(x => x==='inf' ? '∞' : x==='-inf' ? '-∞' : x);

      let lumi_factor = 1;
      if (lumi===0)
        form.querySelector('[name="L"]').value = lumi = r.lumi;
      let data_lumi_text = '';
      if (r.lumi!==lumi) {
        lumi_factor = lumi / r.lumi;
        data_lumi_text = `(scaled from ${r.lumi} ifb)`;
      }
      _id('data_lumi').textContent = data_lumi_text;

      table.childNodes[0].childNodes[7].textContent =
        method==='bkg_rew' ? 'reweighted' :
        method==='sig_reg' ? '[121,129]'  : '';

      for (const [b,bin] of enumerate(r.bins)) {
        const tr = $(table,'tr');
        $(tr,'td').textContent = `[${edges[b]},${edges[b+1]})`;
        let sig_all=0, sig_mid=0;
        let hist = bin.S.hist,
            n3 = hist.length,
            n1 = n3*(121-105)/(160-105),
            n2 = n3*(129-105)/(160-105);
        for (const x of hist)
          sig_all += x;
        for (let i=n1; i<n2; ++i)
          sig_mid += hist[i];
        sig_all *= lumi;
        sig_mid *= lumi;
        const sig =
          method==='bkg_rew' ? sig_all :
          method==='sig_reg' ? sig_mid  : 0;
        $(tr,'td').textContent = sig_all.toFixed(2);
        $(tr,'td').textContent = sig_mid.toFixed(2);
        $(tr,'td').textContent = (100*bin.S.rootw2*lumi/sig).toFixed(2)+'%';
        $(tr,'td').textContent = (100*Math.sqrt(sig_all)/sig).toFixed(2)+'%';

        let Lbkg=0, Rbkg=0;
        hist = bin.B.hist;
        n3 = hist.length*(160-105)/(121-105+160-129);
        n1 = n3*(121-105)/(160-105);
        n2 = n3 - n3*(129-121)/(160-105);
        let i=0;
        for (; i<n1; ++i) Lbkg += hist[i];
        for (; i<n2; ++i) Rbkg += hist[i];
        if (lumi_factor!==1) {
          Lbkg = (Lbkg*lumi_factor).toFixed(2);
          Rbkg = (Rbkg*lumi_factor).toFixed(2);
        }
        $(tr,'td').textContent = Lbkg;
        $(tr,'td').textContent = Rbkg;
        const bkg = bin.B.bkg[ method==='bkg_rew' ? 1 : 0 ];
        $(tr,'td').textContent = bkg.toFixed(2);
        $(tr,'td');
        $(tr,'td').textContent = (100*Math.sqrt(bkg)/bkg).toFixed(2)+'%';
        let signif = sig/Math.sqrt(sig+bkg);
        $(tr,'td',{style:{'font-weight':'bold',color:
          '#006600'
        }}).textContent = signif.toFixed(2);
        signif = Math.sqrt(2*( (sig+bkg)*Math.log(1+sig/bkg)-sig ));
        $(tr,'td',{style:{'font-weight':'bold',color:
          '#006600'
        }}).textContent = signif.toFixed(2);
        $(tr,'td').textContent = (100*sig/(sig+bkg)).toFixed(2)+'%';
        $(tr,'td');
      }
      resize_table();
    });
  }
  form.addEventListener('submit',submit);
  submit();
});

window.addEventListener('resize',resize_table);
