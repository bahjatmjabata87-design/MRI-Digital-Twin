// MRI Digital Twin v1.1 – Positioning & Planning module (brain)
// Educational model only. Head anatomy: MNI ICBM152 2009a GM/WM maps (CSF and scalp derived), 2 mm grid.
(function(){
'use strict';
const VOL={nx:87,ny:105,nz:86,vox:2,origin:[-85.5,-121.5,-73.5]};
const NV=VOL.nx*VOL.ny*VOL.nz;
let D=null; // Float32Array, 4 channels interleaved (gm,wm,csf,fat), 0..1
const $=id=>document.getElementById(id);

// ---------- tissue model ----------
// Brain: same values as the contrast lab. Fat (scalp): 3 T de Bazelaire et al. Radiology 2004 (382/68 ms); 1.5 T approximate.
const TT={'1.5':{gm:[1150,90],wm:[680,70],csf:[4000,2000],fat:[340,75]},
          '3':{gm:[1331,110],wm:[832,80],csf:[4300,2000],fat:[382,68]}};
const TPD={gm:0.80,wm:0.69,csf:1.0,fat:0.90};
// approximate landmarks located on the template (mm, MNI space)
const AC=[0,2,-4], PC=[0,-28,-2];
const BRAIN_LOW=-62; // approximate foramen magnum level on the template (bottom of cerebellum)

// ---------- small vector helpers ----------
const dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
const add=(a,b)=>[a[0]+b[0],a[1]+b[1],a[2]+b[2]];
const sub=(a,b)=>[a[0]-b[0],a[1]-b[1],a[2]-b[2]];
const mul=(a,s)=>[a[0]*s,a[1]*s,a[2]*s];
const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const norm=a=>{const l=Math.hypot(...a)||1;return mul(a,1/l);};
const mv=(M,v)=>[M[0][0]*v[0]+M[0][1]*v[1]+M[0][2]*v[2],M[1][0]*v[0]+M[1][1]*v[1]+M[1][2]*v[2],M[2][0]*v[0]+M[2][1]*v[1]+M[2][2]*v[2]];
const mtv=(M,v)=>[M[0][0]*v[0]+M[1][0]*v[1]+M[2][0]*v[2],M[0][1]*v[0]+M[1][1]*v[1]+M[2][1]*v[2],M[0][2]*v[0]+M[1][2]*v[1]+M[2][2]*v[2]];
const mm=(A,B)=>A.map((r,i)=>[0,1,2].map(j=>A[i][0]*B[0][j]+A[i][1]*B[1][j]+A[i][2]*B[2][j]));
function rot(axis,t){const [x,y,z]=norm(axis),c=Math.cos(t),s=Math.sin(t),C=1-c;
  return [[c+x*x*C,x*y*C-z*s,x*z*C+y*s],[y*x*C+z*s,c+y*y*C,y*z*C-x*s],[z*x*C-y*s,z*y*C+x*s,c+z*z*C]];}
const I3=()=>[[1,0,0],[0,1,0],[0,0,1]];
const deg=r=>r*180/Math.PI, rad=d=>d*Math.PI/180;

// ---------- volume sampling (head coordinates, mm) ----------
const S4=new Float32Array(4);
function sampleHead(q){
  const fx=(q[0]-VOL.origin[0])/2, fy=(q[1]-VOL.origin[1])/2, fz=(q[2]-VOL.origin[2])/2;
  S4[0]=S4[1]=S4[2]=S4[3]=0;
  if(fx<0||fy<0||fz<0||fx>=VOL.nx-1||fy>=VOL.ny-1||fz>=VOL.nz-1) return S4;
  const x=fx|0,y=fy|0,z=fz|0,ax=fx-x,ay=fy-y,az=fz-z;
  const sx=4,sy=4*VOL.nx,sz=4*VOL.nx*VOL.ny, b=x*sx+y*sy+z*sz;
  const w000=(1-ax)*(1-ay)*(1-az),w100=ax*(1-ay)*(1-az),w010=(1-ax)*ay*(1-az),w110=ax*ay*(1-az),
        w001=(1-ax)*(1-ay)*az,w101=ax*(1-ay)*az,w011=(1-ax)*ay*az,w111=ax*ay*az;
  for(let c=0;c<4;c++){
    S4[c]=D[b+c]*w000+D[b+sx+c]*w100+D[b+sy+c]*w010+D[b+sx+sy+c]*w110+
          D[b+sz+c]*w001+D[b+sx+sz+c]*w101+D[b+sy+sz+c]*w011+D[b+sx+sy+sz+c]*w111;
  }
  return S4;
}

// ---------- state ----------
const ST={
  pos:{dir:'hf',pose:'supine',coil:'head',pads:true,ears:true,land:'nasion'},
  moved:false, tilt:I3(), tiltInfo:[0,0],
  L:[0,-15,-30],      // head point placed at isocentre (from laser landmark)
  B0:'1.5',
  task:'ax',
  blk:{type:'ax',R:I3(),c:[0,0,0]},
  p:{TR:4000,TE:100,ETL:16,NEX:1,N:20,thk:5,gap:1,FOV:230,Nf:256,Np:224,phase:1},
  locOff:{sag:1,cor:1,ax:1}, cur:0, scanned:false
};
const LAND={nasion:{z:-30,t:'Nasion'},vertex:{z:94,t:'Vertex (top of the head)'},chin:{z:-110,t:'Chin (mental point)'},sternal:{z:-230,t:'Sternal notch'}};
const TASKS={
  ax:{t:'Routine axial brain (T2 FSE, AC–PC line)',type:'ax',seq:'t2'},
  cor:{t:'Coronal brain (T2 FSE, perpendicular to AC–PC)',type:'cor',seq:'t2'},
  sag:{t:'Sagittal brain (T1 SE, true midline)',type:'sag',seq:'t1'}
};
const SEQ={t1:{TR:500,TE:15,ETL:1,lab:'T1 SE'},t2:{TR:4000,TE:100,ETL:16,lab:'T2 FSE'},pd:{TR:3000,TE:15,ETL:8,lab:'PD FSE'}};

// base frames (scanner = patient coordinates: +x patient left? no: +x = patient RIGHT is MNI convention;
// display is radiological: screen right = patient left = -x)
const BASE={
  ax:{n:[0,0,1],e1:[-1,0,0],e2:[0,-1,0],lab:['R','L','A','P']},
  cor:{n:[0,1,0],e1:[-1,0,0],e2:[0,0,-1],lab:['R','L','S','I']},
  sag:{n:[1,0,0],e1:[0,-1,0],e2:[0,0,-1],lab:['A','P','S','I']}
};
function frame(){const b=BASE[ST.blk.type],R=ST.blk.R;return {n:mv(R,b.n),e1:mv(R,b.e1),e2:mv(R,b.e2)};}
const toHead=p=>add(mtv(ST.tilt,p),ST.L);           // scanner -> head
const toScan=q=>mv(ST.tilt,sub(q,ST.L));            // head -> scanner

// ---------- signal ----------
function tissueSignals(){
  const T=TT[ST.B0],p=ST.p,o={};
  for(const k of ['gm','wm','csf','fat']){
    let [T1,T2]=T[k]; if(k==='fat'&&p.ETL>1) T2*=2.5; // J-coupling: fat stays bright on fast spin echo
    o[k]=TPD[k]*(1-Math.exp(-p.TR/T1))*Math.exp(-p.TE/T2);
  }
  return o;
}
function weighting(TR,TE){
  if(TR<=800&&TE<=30) return 'T1-weighted';
  if(TR>=2000&&TE>=80) return 'T2-weighted';
  if(TR>=2000&&TE<=30) return 'PD-weighted';
  if(TR<=800&&TE>=80) return 'Poor combination (low signal)';
  return 'Mixed weighting';
}
const REF={vox:(230/256)*(230/256)*5,Np:256,NEX:1};
function physics(){
  const p=ST.p, px=p.FOV/p.Nf, py=p.FOV/p.Np;
  const vox=px*py*p.thk;
  const coil=ST.pos.coil==='head'?1:0.35;
  const rel=(vox/REF.vox)*Math.sqrt(p.Np*p.NEX/(REF.Np*REF.NEX))*coil*(ST.B0==='3'?2:1);
  const ESP=10, perSlice=Math.max(p.TE,p.ETL*ESP)+10;
  const maxPerTR=Math.max(1,Math.floor(p.TR/perSlice));
  const conc=Math.ceil(p.N/maxPerTR);
  const tsec=p.TR*p.Np*p.NEX/p.ETL*conc/1000;
  const cov=p.N*p.thk+(p.N-1)*p.gap;
  return {px,py,vox,rel,maxPerTR,conc,tsec,cov};
}
const fmtT=s=>{s=Math.round(s);return Math.floor(s/60)+' min '+String(s%60).padStart(2,'0')+' s';};

// ---------- localizer views ----------
const LV={
  sag:{u:[0,-1,0],w:[0,0,-1],m:[1,0,0],lab:['A','P','S','I'],t:'Sagittal localizer'},
  cor:{u:[-1,0,0],w:[0,0,-1],m:[0,1,0],lab:['R','L','S','I'],t:'Coronal localizer'},
  ax:{u:[-1,0,0],w:[0,-1,0],m:[0,0,1],lab:['R','L','A','P'],t:'Axial localizer'}
};
const LFOV=280, LPX=300, LOFFS=[-15,0,15];
const locCache={};
function locCenter(k){return mul(LV[k].m,LOFFS[ST.locOff[k]]);}
function renderLocalizer(k){
  const v=LV[k],o=locCenter(k),s=LFOV/LPX,img=new ImageData(LPX,LPX);
  for(let j=0;j<LPX;j++)for(let i=0;i<LPX;i++){
    const a=(i-LPX/2+.5)*s,b=(j-LPX/2+.5)*s;
    const p=[o[0]+a*v.u[0]+b*v.w[0],o[1]+a*v.u[1]+b*v.w[1],o[2]+a*v.u[2]+b*v.w[2]];
    const t=sampleHead(toHead(p));
    // fast gradient-echo style localizer contrast
    let val=(t[0]*0.55+t[1]*0.8+t[2]*0.18+t[3]*0.95)*250+Math.random()*8;
    if(ST.pos.coil!=='head') val*=0.75;
    const c=Math.max(0,Math.min(255,val)),q=(j*LPX+i)*4;
    img.data[q]=img.data[q+1]=img.data[q+2]=c; img.data[q+3]=255;
  }
  locCache[k]=img;
}
function toScreen(k,p){const v=LV[k],o=locCenter(k),d=sub(p,o),s=LPX/LFOV;return [LPX/2+dot(d,v.u)*s,LPX/2+dot(d,v.w)*s];}
function fromScreen(k,x,y){const v=LV[k],o=locCenter(k),s=LFOV/LPX;return add(o,add(mul(v.u,(x-LPX/2)*s),mul(v.w,(y-LPX/2)*s)));}

function sliceCenters(){const f=frame(),p=ST.p,out=[];for(let i=0;i<p.N;i++)out.push(add(ST.blk.c,mul(f.n,(i-(p.N-1)/2)*(p.thk+p.gap))));return out;}
function segInView(k,ck,f){
  const v=LV[k],m=v.m,o=locCenter(k),n=f.n,F=ST.p.FOV/2;
  let d=cross(n,m); const dl=Math.hypot(...d); if(dl<0.2) return null;
  const p0=mul(add(mul(cross(m,d),dot(n,ck)),mul(cross(d,n),dot(m,o))),1/(dl*dl));
  d=mul(d,1/dl); let t0=-1e4,t1=1e4;
  for(const e of [f.e1,f.e2]){const a=dot(sub(p0,ck),e),b=dot(d,e);
    if(Math.abs(b)<1e-6){if(Math.abs(a)>F) return null; continue;}
    let lo=(-F-a)/b,hi=(F-a)/b; if(lo>hi)[lo,hi]=[hi,lo]; t0=Math.max(t0,lo); t1=Math.min(t1,hi);}
  if(t0>=t1) return null;
  return [add(p0,mul(d,t0)),add(p0,mul(d,t1))];
}
function drawLocalizer(k){
  const cv=$('loc-'+k),g=cv.getContext('2d');
  g.fillStyle='#000';g.fillRect(0,0,LPX,LPX);
  if(!ST.moved||!locCache[k]){g.fillStyle='#8a94a6';g.font='14px system-ui';g.textAlign='center';g.fillText('Position the patient first',LPX/2,LPX/2);g.textAlign='left';return;}
  g.putImageData(locCache[k],0,0);
  const v=LV[k];
  g.fillStyle='rgba(255,255,255,.85)';g.font='bold 12px system-ui';
  g.fillText(v.lab[0],6,LPX/2+4);g.fillText(v.lab[1],LPX-14,LPX/2+4);g.fillText(v.lab[2],LPX/2-4,14);g.fillText(v.lab[3],LPX/2-3,LPX-6);
  g.font='11px system-ui';g.fillText('Slice '+(ST.locOff[k]+1)+'/3 · '+(LOFFS[ST.locOff[k]]>=0?'+':'')+LOFFS[ST.locOff[k]]+' mm',6,LPX-6);
  if($('showacpc').checked&&k==='sag'){
    const a=toScreen(k,toScan(AC)),b=toScreen(k,toScan(PC)),dx=a[0]-b[0],dy=a[1]-b[1],L=Math.hypot(dx,dy);
    g.strokeStyle='#FFD166';g.setLineDash([6,4]);g.lineWidth=1.5;g.beginPath();g.moveTo(b[0]-dx/L*120,b[1]-dy/L*120);g.lineTo(a[0]+dx/L*120,a[1]+dy/L*120);g.stroke();g.setLineDash([]);
    for(const [p,t] of [[a,'AC'],[b,'PC']]){g.fillStyle='#FFD166';g.beginPath();g.arc(p[0],p[1],3.5,0,7);g.fill();g.font='bold 11px system-ui';g.fillText(t,p[0]-8,p[1]-7);}
  }
  const f=frame(),cs=sliceCenters();
  const segs=cs.map(c=>segInView(k,c,f));
  if(segs.every(s=>!s)){ // view parallel to the slices: show FOV square of the central slice
    const F=ST.p.FOV/2,c=ST.blk.c,pts=[[-1,-1],[1,-1],[1,1],[-1,1]].map(([a,b])=>toScreen(k,add(c,add(mul(f.e1,a*F),mul(f.e2,b*F)))));
    g.strokeStyle='#00B4D8';g.lineWidth=2;g.setLineDash([5,4]);g.beginPath();pts.forEach((p,i)=>i?g.lineTo(...p):g.moveTo(...p));g.closePath();g.stroke();g.setLineDash([]);
  } else {
    const first=segs.find(s=>s),last=[...segs].reverse().find(s=>s);
    if(first&&last){const P=[first[0],first[1],last[1],last[0]].map(p=>toScreen(k,p));
      g.fillStyle='rgba(0,180,216,.12)';g.beginPath();P.forEach((p,i)=>i?g.lineTo(...p):g.moveTo(...p));g.closePath();g.fill();}
    segs.forEach((s,i)=>{if(!s)return;const a=toScreen(k,s[0]),b=toScreen(k,s[1]);
      const curr=ST.scanned&&i===ST.cur;
      g.strokeStyle=curr?'#FFD166':'#00B4D8';g.lineWidth=curr?2.5:1;g.beginPath();g.moveTo(...a);g.lineTo(...b);g.stroke();});
  }
  // centre + rotation handle
  const c=toScreen(k,ST.blk.c),h=handlePos(k);
  g.fillStyle='#EF476F';g.beginPath();g.arc(c[0],c[1],4,0,7);g.fill();
  if(h){g.strokeStyle='#EF476F';g.lineWidth=1.5;g.beginPath();g.moveTo(...c);g.lineTo(...h);g.stroke();
    g.fillStyle='#fff';g.beginPath();g.arc(h[0],h[1],7,0,7);g.fill();g.strokeStyle='#EF476F';g.lineWidth=2;g.stroke();
    g.fillStyle='#EF476F';g.font='bold 10px system-ui';g.fillText('⟳',h[0]-4,h[1]+4);}
}
function handlePos(k){
  const f=frame(),s=segInView(k,ST.blk.c,f);
  if(!s) return null;
  const a=toScreen(k,s[0]),b=toScreen(k,s[1]),c=toScreen(k,ST.blk.c);
  let dx=b[0]-a[0],dy=b[1]-a[1];const L=Math.hypot(dx,dy)||1;dx/=L;dy/=L;
  if(dx<0||(Math.abs(dx)<1e-3&&dy<0)){dx=-dx;dy=-dy;}
  const r=Math.min(L/2+14,135);
  return [c[0]+dx*r,c[1]+dy*r];
}
function drawAllLoc(){['sag','cor','ax'].forEach(drawLocalizer);}

// ---------- pointer interaction on localizers ----------
function bindLocalizer(k){
  const cv=$('loc-'+k);let mode=null,last=null;
  const pt=e=>{const r=cv.getBoundingClientRect();return [(e.clientX-r.left)*LPX/r.width,(e.clientY-r.top)*LPX/r.height];};
  cv.addEventListener('pointerdown',e=>{
    if(!ST.moved) return; const p=pt(e),h=handlePos(k);
    mode=(h&&Math.hypot(p[0]-h[0],p[1]-h[1])<14)?'rot':'move'; last=p; cv.setPointerCapture(e.pointerId); e.preventDefault();
  });
  cv.addEventListener('pointermove',e=>{
    if(!mode) return; const p=pt(e);
    if(mode==='move'){const a=fromScreen(k,...last),b=fromScreen(k,...p);ST.blk.c=add(ST.blk.c,sub(b,a));}
    else{const c=toScreen(k,ST.blk.c);
      const a0=Math.atan2(last[1]-c[1],last[0]-c[0]),a1=Math.atan2(p[1]-c[1],p[0]-c[0]);
      // screen angle -> rotation about the view normal (u x w)
      const ax=cross(LV[k].u,LV[k].w);ST.blk.R=mm(rot(ax,a1-a0),ST.blk.R);}
    last=p; changed();
  });
  const up=()=>{mode=null;};
  cv.addEventListener('pointerup',up);cv.addEventListener('pointercancel',up);
}

// ---------- result images ----------
const RPX=320;
function txt(g,t,x,y){g.save();g.shadowColor='#000';g.shadowBlur=4;g.lineWidth=3;g.strokeStyle='rgba(0,0,0,.85)';g.strokeText(t,x,y);g.fillText(t,x,y);g.restore();}
function renderResult(){
  const cv=$('res'),g=cv.getContext('2d');
  g.fillStyle='#000';g.fillRect(0,0,RPX,RPX);
  if(!ST.scanned){g.fillStyle='#8a94a6';g.font='14px system-ui';g.textAlign='center';g.fillText('Plan the slices, then press Scan',RPX/2,RPX/2);g.textAlign='left';return;}
  const p=ST.p,f=frame(),cs=sliceCenters(),ck=cs[Math.min(ST.cur,cs.length-1)];
  const sig=tissueSignals(),ph=physics();
  const sd=0.012/ph.rel, K=255/0.45;
  const Nf=p.Nf,Np=p.Np, phaseIsE1=p.phase===1;
  const nx=phaseIsE1?Np:Nf, ny=phaseIsE1?Nf:Np; // e1 = image columns, e2 = rows
  const F=p.FOV, offs=[-p.thk/3,0,p.thk/3];
  const raw=new Float32Array(nx*ny);
  const wrapDir=phaseIsE1?f.e1:f.e2;
  for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){
    const a=((i+.5)/nx-.5)*F,b=((j+.5)/ny-.5)*F;let v=0;
    for(const w of [-1,0,1]){ // aliasing: tissue one FOV away along the phase direction folds in
      const base=add(ck,add(add(mul(f.e1,a),mul(f.e2,b)),mul(wrapDir,w*F)));
      for(const o of offs){const ps=add(base,mul(f.n,o));const t=sampleHead(toHead(ps));
        let s=t[0]*sig.gm+t[1]*sig.wm+t[2]*sig.csf+t[3]*sig.fat;
        const r=Math.hypot(...ps); s*=Math.exp(-Math.pow(r/210,6)); // field homogeneity / coil falloff far from isocentre
        v+=s/3;}
    }
    raw[j*nx+i]=v;
  }
  let out=raw;
  if(!ST.pos.pads){ // head motion -> ghosts along the phase-encoding direction
    out=new Float32Array(raw.length);const sh=Math.round((phaseIsE1?nx:ny)/5);
    for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){let v=raw[j*nx+i];
      for(const s of [-sh,sh]){if(phaseIsE1){const ii=(i+s+nx)%nx;v+=0.28*raw[j*nx+ii];}else{const jj=(j+s+ny)%ny;v+=0.28*raw[jj*nx+i];}}
      out[j*nx+i]=v;}
  }
  const off=document.createElement('canvas');off.width=nx;off.height=ny;const og=off.getContext('2d'),img=og.createImageData(nx,ny);
  for(let q=0;q<nx*ny;q++){const v=Math.hypot(out[q]+gauss()*sd,gauss()*sd)*K,c=Math.max(0,Math.min(255,v));
    img.data[q*4]=img.data[q*4+1]=img.data[q*4+2]=c;img.data[q*4+3]=255;}
  og.putImageData(img,0,0);
  g.imageSmoothingEnabled=false;g.drawImage(off,0,0,RPX,RPX);
  const lab=BASE[ST.blk.type].lab;
  g.fillStyle='rgba(255,255,255,.85)';g.font='bold 12px system-ui';
  g.fillText(lab[0],6,RPX/2+4);g.fillText(lab[1],RPX-14,RPX/2+4);g.fillText(lab[2],RPX/2-4,14);g.fillText(lab[3],RPX/2-3,RPX-6);
  g.font='11px system-ui';
  txt(g,SEQLAB()+' · '+ST.B0+' T · TR '+p.TR+' · TE '+p.TE,8,30);
  txt(g,'Slice '+(ST.cur+1)+' / '+p.N+' · '+p.thk+' mm · matrix '+p.Nf+'×'+p.Np,8,RPX-22);
}
function gauss(){let u=0,v=0;while(!u)u=Math.random();while(!v)v=Math.random();return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v);}
const SEQLAB=()=>ST.p.ETL>1?'FSE (ETL '+ST.p.ETL+')':'SE';

// ---------- task checks ----------
let BRAINPTS=null, SCALPPTS=null;
function buildPoints(){
  BRAINPTS=[];SCALPPTS=[];
  for(let z=0;z<VOL.nz;z+=2)for(let y=0;y<VOL.ny;y+=2)for(let x=0;x<VOL.nx;x+=2){
    const b=4*(x+VOL.nx*(y+VOL.ny*z)),q=[VOL.origin[0]+2*x,VOL.origin[1]+2*y,VOL.origin[2]+2*z];
    if(D[b]+D[b+1]>0.5&&q[2]>=BRAIN_LOW) BRAINPTS.push(q);
    if(D[b+3]>0.5) SCALPPTS.push(q);
  }
}
function inBlock(p,f,half,F){const d=sub(p,ST.blk.c);return Math.abs(dot(d,f.n))<=half&&Math.abs(dot(d,f.e1))<=F&&Math.abs(dot(d,f.e2))<=F;}
function angleErr(){
  const nh=norm(mtv(ST.tilt,frame().n)); // block normal in head coordinates
  const l=norm(sub(PC,AC)), xh=[1,0,0];
  let ref; if(ST.blk.type==='ax') ref=norm(cross(xh,l)); else if(ST.blk.type==='cor') ref=l; else ref=xh;
  if(ref[2]<0&&ST.blk.type==='ax') ref=mul(ref,-1);
  if(dot(nh,ref)<0) ref=mul(ref,-1);
  // components: tilt seen on the sagittal (y-z) and coronal/axial planes
  const pitch=deg(Math.atan2(nh[1]*ref[2]-nh[2]*ref[1], nh[1]*ref[1]+nh[2]*ref[2]));
  const roll=deg(Math.atan2(nh[0]*ref[2]-nh[2]*ref[0], nh[0]*ref[0]+nh[2]*ref[2]));
  const yaw=deg(Math.atan2(nh[0]*ref[1]-nh[1]*ref[0], nh[0]*ref[0]+nh[1]*ref[1]));
  return {total:deg(Math.acos(Math.min(1,Math.abs(dot(nh,ref))))),pitch,roll,yaw};
}
function runChecks(){
  const out=[],P=ST.pos,p=ST.p,ph=physics(),f=frame(),T=TASKS[ST.task];
  const ok=(c,t,why)=>out.push({c,t,why});
  ok(P.dir==='hf'&&P.pose==='supine','Head first, supine','Standard for brain: comfortable, airway visible, and the head reaches the isocentre first.');
  ok(P.coil==='head','Dedicated head coil','A close-fitting receive coil gives much higher SNR than the body coil.');
  ok(P.pads,'Foam pads for immobilisation','Pads reduce head motion, which causes ghosting along the phase-encoding direction.');
  ok(P.ears,'Hearing protection','Gradient noise is loud; earplugs or headphones are required.');
  ok(P.land==='nasion','Laser landmark on the nasion','Centres the head at the isocentre, where the field is most uniform.');
  ok(ST.blk.type===T.type,'Orientation: '+({ax:'axial',cor:'coronal',sag:'sagittal'})[T.type],'The task asks for this plane.');
  const ae=angleErr();
  const angTxt={ax:'Axial slices parallel to the AC–PC line',cor:'Coronal slices perpendicular to the AC–PC line',sag:'Sagittal slices parallel to the midline (falx)'}[ST.blk.type];
  const cand={ax:[['sagittal',ae.pitch],['coronal',ae.roll]],cor:[['sagittal',ae.pitch],['axial',ae.yaw]],sag:[['coronal',ae.roll],['axial',ae.yaw]]}[ST.blk.type].sort((a,b)=>Math.abs(b[1])-Math.abs(a[1]));
  ok(ae.total<=5,angTxt+' (error '+ae.total.toFixed(1)+'°)','Within 5° gives comparable, symmetrical images between visits.'+(ae.total>5?' Correct the angle mainly on the '+cand[0][0]+' localizer.':''));
  const half=ph.cov/2,F=p.FOV/2;let inn=0,top=0,bot=0;
  for(const q of BRAINPTS){const s=toScan(q);if(inBlock(s,f,half,F))inn++;else{const d=dot(sub(s,ST.blk.c),f.n);if(d>0)top++;else bot++;}}
  const pct=100*inn/BRAINPTS.length;
  ok(pct>=98,'Whole brain covered ('+pct.toFixed(1)+'%)',ST.blk.type==='ax'?'Cover from the vertex to the foramen magnum.'+(pct<98?(top>bot?' Some brain is missed at the top (vertex) or one edge – move up or add slices.':' Some brain is missed at the bottom or one edge – move down or add slices.'):''):'Cover the whole brain from side to side or front to back, inside the FOV.');
  let wrap=0;const ax=p.phase===1?f.e1:f.e2;
  for(const q of SCALPPTS){const s=toScan(q),d=sub(s,ST.blk.c);if(Math.abs(dot(d,f.n))<=half&&Math.abs(dot(d,ax))>F)wrap++;}
  ok(wrap<3,'No wrap-around (aliasing)','The head must fit inside the FOV along the phase-encoding direction, or tissue outside folds into the image.');
  const want=T.seq==='t2'?'T2-weighted':'T1-weighted';
  ok(weighting(p.TR,p.TE)===want,'Contrast: '+want,'Task asks for '+want+'. T1: short TR/short TE; T2: long TR/long TE.');
  ok(p.thk>=3&&p.thk<=6,'Slice thickness 3–6 mm','Thicker slices raise SNR but blur small lesions (partial volume).');
  ok(p.gap<=p.thk*0.3,'Gap ≤ 30% of thickness','Large gaps can miss small lesions.');
  ok(ph.px<=1.0&&ph.py<=1.2,'In-plane pixel ≤ 1.0 × 1.2 mm ('+ph.px.toFixed(2)+' × '+ph.py.toFixed(2)+')','Enough resolution for brain anatomy.');
  ok(ph.rel>=0.8,'Relative SNR ≥ 80% ('+Math.round(ph.rel*100)+'%)','Smaller voxels or fewer averages lower SNR.');
  ok(ph.tsec<=360,'Scan time ≤ 6 min ('+fmtT(ph.tsec)+')','Long scans increase the chance of patient motion.');
  return out;
}
function showChecks(){
  const r=runChecks(),n=r.filter(x=>x.c).length;
  $('chk').innerHTML='<div class="pill">Score '+n+' / '+r.length+'</div>'+r.map(x=>`<div class="ck ${x.c?'y':'n'}"><b>${x.c?'✓':'✗'} ${x.t}</b><span>${x.why}</span></div>`).join('');
}

// ---------- positioning diagram ----------
let anim=0;
function drawPos(){
  const cv=$('posc'),g=cv.getContext('2d'),W=cv.width,H=cv.height,P=ST.pos;
  g.clearRect(0,0,W,H);
  // bore
  const bx=W*0.42,bw=W*0.30;
  g.fillStyle='#dfe6f2';g.fillRect(0,0,W,H);
  g.fillStyle='#b9c4d8';g.beginPath();g.roundRect(bx,20,bw,H-60,26);g.fill();
  g.fillStyle='#eef2f9';g.beginPath();g.roundRect(bx+10,62,bw-20,H-150,40);g.fill();
  g.fillStyle='#5A6472';g.font='12px system-ui';g.fillText('Magnet bore',bx+bw/2-34,40);
  const iso=bx+bw/2;
  g.strokeStyle='#EF476F';g.setLineDash([4,4]);g.beginPath();g.moveTo(iso,62);g.lineTo(iso,H-88);g.stroke();g.setLineDash([]);
  g.fillStyle='#EF476F';g.fillText('Isocentre',iso-24,H-72);
  // table + patient (slides in with anim 0..1)
  const lz=LAND[P.land].z; // mm, head coordinate placed at laser
  const scale=0.9, travel=(iso-(W*0.20))*anim;
  const headTop=H-110, tableY=H-92;
  g.fillStyle='#8a94a6';g.fillRect(20,tableY,W-40,10);
  const hf=P.dir==='hf', sgn=hf?1:-1;
  // head x (centre of brain) before travel: laser line at W*0.20 marks the landmark
  const laserX=W*0.20;
  const headCx=laserX+travel+sgn*(10-lz)*scale*0.55; // brain centre about +10 mm (MNI z)
  g.save();
  // body
  const bodyX0=hf?headCx-sgn*0:headCx;
  g.fillStyle='#F5C6A5';
  const bodyStart=headCx-sgn*(58), bodyEnd=bodyStart-sgn*340;
  g.beginPath();g.roundRect(Math.min(bodyStart,bodyEnd),headTop+18,Math.abs(bodyEnd-bodyStart),40,16);g.fill();
  // head
  g.beginPath();g.ellipse(headCx,headTop+30,46,32,0,0,7);g.fill();
  // nose / face direction
  const faceUp=P.pose==='supine';
  g.beginPath();g.moveTo(headCx+sgn*(-6),headTop+30+(faceUp?-30:30));g.lineTo(headCx+sgn*10,headTop+30+(faceUp?-44:44));g.lineTo(headCx+sgn*18,headTop+30+(faceUp?-28:28));g.fill();
  g.fillStyle='#1B1F3B';g.beginPath();g.arc(headCx+sgn*6,headTop+30+(faceUp?-18:18),3,0,7);g.fill();
  if(P.pads){g.fillStyle='#06D6A0';g.fillRect(headCx-40,headTop-6,80,7);g.fillRect(headCx-40,headTop+60,80,7);}
  if(P.ears){g.fillStyle='#FFD166';g.beginPath();g.arc(headCx-sgn*8,headTop+34,6,0,7);g.fill();}
  if(P.coil==='head'){g.strokeStyle='#3A0CA3';g.lineWidth=4;g.beginPath();g.roundRect(headCx-58,headTop-22,116,104,18);g.stroke();
    g.lineWidth=2;for(let i=-40;i<=40;i+=20){g.beginPath();g.moveTo(headCx+i,headTop-22);g.lineTo(headCx+i,headTop+82);g.stroke();}}
  g.restore();
  // laser (fixed on the gantry front, before travel)
  if(anim<0.02){g.strokeStyle='#ff2d2d';g.lineWidth=2;g.beginPath();g.moveTo(laserX,8);g.lineTo(laserX,tableY);g.stroke();
    g.fillStyle='#ff2d2d';g.font='bold 12px system-ui';g.fillText('Laser: '+LAND[P.land].t,laserX+6,16);}
  g.fillStyle='#1B1F3B';g.font='12px system-ui';
  g.fillText((hf?'Head first':'Feet first')+' · '+(P.pose==='supine'?'Supine':'Prone')+' · '+(P.coil==='head'?'Head coil':'Body coil'),12,H-12);
}
function moveIn(){
  ST.moved=false;anim=0;const t0=performance.now();$('movebtn').disabled=true;
  const step=now=>{anim=Math.min(1,(now-t0)/1400);drawPos();if(anim<1)requestAnimationFrame(step);else{
    ST.L=[0,-15,LAND[ST.pos.land].z];
    ST.moved=true;['sag','cor','ax'].forEach(renderLocalizer);$('movebtn').disabled=false;
    $('posmsg').innerHTML='Patient at the isocentre. Localizers acquired – plan your slices in Step 2.';changed();}};
  requestAnimationFrame(step);
}
function newPatient(){
  const pitch=(Math.random()<.5?-1:1)*(6+Math.random()*8), roll=(Math.random()<.5?-1:1)*(1+Math.random()*3), yaw=(Math.random()<.5?-1:1)*(Math.random()*3);
  ST.tilt=mm(rot([0,0,1],rad(yaw)),mm(rot([0,1,0],rad(roll)),rot([1,0,0],rad(pitch))));
  ST.tiltInfo=[pitch,roll];
  resetPlan(); ST.moved=false; ST.scanned=false; anim=0; drawPos(); changed();
  $('posmsg').innerHTML='New patient on the table. Like real patients, the head is not perfectly straight – you will correct this with the slice angle.';
}
function resetPlan(){const T=TASKS[ST.task];ST.blk={type:T.type,R:I3(),c:[0,0,0]};ST.cur=0;}
function autoPlan(){
  // demo: align to the reference line, centre on the brain, apply task protocol
  const T=TASKS[ST.task];ST.blk.type=T.type;
  const l=norm(sub(PC,AC)),xh=[1,0,0];
  let nH=T.type==='ax'?norm(cross(xh,l)):T.type==='cor'?l:xh; if(T.type==='ax'&&nH[2]<0) nH=mul(nH,-1); if(T.type==='cor'&&nH[1]<0)nH=mul(nH,-1);
  const e1H=T.type==='sag'?[0,-1,0]:[-1,0,0];
  // build head-frame basis then map to scanner
  const b=BASE[T.type], hand=Math.sign(dot(cross(b.n,b.e1),b.e2));
  const nS=norm(mv(ST.tilt,nH)); let e1S=mv(ST.tilt,e1H); e1S=norm(sub(e1S,mul(nS,dot(e1S,nS)))); const e2S=mul(cross(nS,e1S),hand); // R maps base frame to (nS,e1S,e2S)
  const M=[nS,e1S,e2S], Bm=[b.n,b.e1,b.e2];
  const R=[[0,0,0],[0,0,0],[0,0,0]];
  for(let i=0;i<3;i++)for(let j=0;j<3;j++){let s=0;for(let k=0;k<3;k++)s+=M[k][i]*Bm[k][j];R[i][j]=s;}
  ST.blk.R=R;
  const ctrH=T.type==='ax'?[0,-18,(82+BRAIN_LOW)/2]:T.type==='cor'?[0,-17,8]:[0,-17,8];
  ST.blk.c=toScan(ctrH);
  const pr=T.seq==='t2'?SEQ.t2:SEQ.t1;Object.assign(ST.p,{TR:pr.TR,TE:pr.TE,ETL:pr.ETL,NEX:T.type==='sag'?2:1,FOV:230,Nf:256,Np:224});
  if(T.type==='ax'){ST.p.N=24;ST.p.thk=5;ST.p.gap=1;ST.p.phase=2;}
  else if(T.type==='cor'){ST.p.N=26;ST.p.thk=5;ST.p.gap=1;ST.p.phase=1;}
  else {ST.p.N=25;ST.p.thk=5;ST.p.gap=1;ST.p.phase=1;}
  syncInputs();changed();
}

// ---------- UI wiring ----------
const PARAMS=[['TR','ms'],['TE','ms'],['ETL',''],['NEX',''],['N',''],['thk','mm'],['gap','mm'],['FOV','mm']];
function syncInputs(){
  for(const [k] of PARAMS){$('pp-'+k).value=ST.p[k];}
  $('pp-Nf').value=ST.p.Nf;$('pp-Np').value=ST.p.Np;$('pp-phase').value=ST.p.phase;
  $('pl-type').value=ST.blk.type;$('pl-task').value=ST.task;
  document.querySelectorAll('[data-pb]').forEach(b=>b.classList.toggle('alt',b.dataset.pb!==ST.B0));
  const P=ST.pos;$('ps-dir').value=P.dir;$('ps-pose').value=P.pose;$('ps-coil').value=P.coil;$('ps-pads').checked=P.pads;$('ps-ears').checked=P.ears;$('ps-land').value=P.land;
  updatePhaseLabels();
}
function updatePhaseLabels(){
  const L={ax:['Right–left','Anterior–posterior'],cor:['Right–left','Superior–inferior'],sag:['Anterior–posterior','Superior–inferior']}[ST.blk.type];
  const s=$('pp-phase');s.options[0].text=L[0];s.options[1].text=L[1];
}
let pend=0;
function changed(){
  const p=ST.p,ph=physics();
  for(const [k,u] of PARAMS) $('pv-'+k).textContent=p[k]+(u?' '+u:'');
  if(ST.cur>=p.N) ST.cur=p.N-1;
  $('rs').max=p.N-1;$('rs').value=ST.cur;
  const ae=angleErr();
  $('readout').innerHTML=`<table>
<tr><td>Contrast</td><td><b>${weighting(p.TR,p.TE)}</b> (${SEQLAB()})</td></tr>
<tr><td>Coverage</td><td><b>${ph.cov.toFixed(0)} mm</b> (${p.N} × ${p.thk} mm + ${p.N-1} × ${p.gap} mm gap)</td></tr>
<tr><td>Voxel size</td><td><b>${ph.px.toFixed(2)} × ${ph.py.toFixed(2)} × ${p.thk} mm</b> (FOV ÷ matrix)</td></tr>
<tr><td>Relative SNR</td><td><b>${Math.round(ph.rel*100)}%</b> <span class="muted">(100% = 230 mm, 256×256, 5 mm, 1 NEX, 1.5 T, head coil)</span></td></tr>
<tr><td>Slices per TR</td><td><b>${Math.min(ph.maxPerTR,p.N)}</b>${ph.conc>1?` – needs <b>${ph.conc}</b> acquisitions (time ×${ph.conc})`:''}</td></tr>
<tr><td>Scan time</td><td><b>${fmtT(ph.tsec)}</b> <span class="muted">= TR × phase steps × NEX ÷ ETL${ph.conc>1?' × acquisitions':''}</span></td></tr>
<tr><td>Block angle</td><td><b>${ST.moved?ae.total.toFixed(1)+'°':'–'}</b> from the ${ST.blk.type==='sag'?'midline':'AC–PC reference'}</td></tr></table>`;
  drawAllLoc();
  if(ST.scanned){clearTimeout(pend);pend=setTimeout(()=>{renderResult();showChecks();},60);} else renderResult();
}
function init(){
  const html=Object.entries(TASKS).map(([k,v])=>`<option value="${k}">${v.t}</option>`).join('');$('pl-task').innerHTML=html;
  for(const [k] of PARAMS){$('pp-'+k).addEventListener('input',e=>{ST.p[k]=+e.target.value;changed();});}
  for(const k of ['Nf','Np','phase']) $('pp-'+k).addEventListener('change',e=>{ST.p[k]=+e.target.value;changed();});
  document.querySelectorAll('[data-pb]').forEach(b=>b.onclick=()=>{ST.B0=b.dataset.pb;syncInputs();changed();});
  document.querySelectorAll('[data-ps]').forEach(b=>b.onclick=()=>{const s=SEQ[b.dataset.ps];Object.assign(ST.p,{TR:s.TR,TE:s.TE,ETL:s.ETL});syncInputs();changed();});
  $('pl-type').onchange=e=>{ST.blk.type=e.target.value;ST.blk.R=I3();updatePhaseLabels();changed();};
  $('pl-task').onchange=e=>{ST.task=e.target.value;resetPlan();syncInputs();ST.scanned=false;$('chk').innerHTML='';changed();};
  const posIn=[['ps-dir','dir'],['ps-pose','pose'],['ps-coil','coil'],['ps-land','land']];
  for(const [id,k] of posIn) $(id).onchange=e=>{ST.pos[k]=e.target.value;ST.moved=false;anim=0;drawPos();$('posmsg').textContent='Setup changed – move the table to the isocentre again.';changed();};
  $('ps-pads').onchange=e=>{ST.pos.pads=e.target.checked;drawPos();changed();};
  $('ps-ears').onchange=e=>{ST.pos.ears=e.target.checked;drawPos();changed();};
  $('movebtn').onclick=moveIn; $('newpt').onclick=newPatient;
  $('resetplan').onclick=()=>{resetPlan();changed();}; $('autoplan').onclick=()=>{if(!ST.moved){$('posmsg').textContent='Move the patient to the isocentre first.';return;}autoPlan();};
  $('showacpc').onchange=drawAllLoc;
  document.querySelectorAll('[data-lo]').forEach(b=>b.onclick=()=>{const [k,d]=b.dataset.lo.split(':');ST.locOff[k]=Math.max(0,Math.min(2,ST.locOff[k]+ +d));drawLocalizer(k);});
  $('scanbtn').onclick=()=>{if(!ST.moved){$('posmsg').textContent='Move the patient to the isocentre first.';return;}ST.scanned=true;ST.cur=Math.floor(ST.p.N/2);changed();$('res').scrollIntoView({behavior:'smooth',block:'center'});};
  $('rs').oninput=e=>{ST.cur=+e.target.value;changed();};
  ['sag','cor','ax'].forEach(bindLocalizer);
  newPatient(); syncInputs(); changed();
}
Promise.resolve().then(()=>{if(!window.BRAIN3D)throw new Error('no data');
  const bin=atob(window.BRAIN3D),u=new Uint8Array(bin.length);for(let i=0;i<bin.length;i++)u[i]=bin.charCodeAt(i);D=new Float32Array(NV*4);
  for(let i=0;i<NV;i++){const a=u[i],b=u[NV+i];D[i*4]=(a&15)/15;D[i*4+1]=(a>>4)/15;D[i*4+2]=(b&15)/15;D[i*4+3]=(b>>4)/15;}
  buildPoints();init();$('planload').style.display='none';
}).catch(e=>{$('planload').textContent='Could not load the 3D head data (brain3d.js). Keep brain3d.js in the same folder as index.html.';});
window.__plan=ST; window.__planHandle=k=>handlePos(k); window.__planErr=()=>angleErr(); window.__planAuto=()=>autoPlan(); window.__planMove=()=>moveIn();
})();
