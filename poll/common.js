// Shared helpers for the participant and presenter pages.
(function () {
  const params = new URLSearchParams(location.search);
  const cfg = window.LIVE_POLL || {};
  // ?sb=http://localhost:8787 points the page at a local stand-in (dev/server.js).
  const BASE = (params.get('sb') || cfg.supabaseUrl || '').replace(/\/+$/, '');
  const KEY = cfg.publishableKey || '';

  // Errors that retrying won't fix.
  const FATAL = /No questions|Missing participant|must be an object/;

  const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

  // Calls one of the database functions in supabase/schema.sql. Network blips
  // and 5xx/429 responses are retried with jittered backoff. Resubmitting is
  // harmless: the newest answer from a participant wins.
  async function rpc(fn, args, tries = 4) {
    const headers = { 'Content-Type': 'application/json', apikey: KEY };
    // Legacy JWT anon keys also go in Authorization; sb_publishable_ keys don't.
    if (KEY.startsWith('eyJ')) headers.Authorization = 'Bearer ' + KEY;
    let lastErr;
    for (let i = 0; i < tries; i++) {
      try {
        const res = await fetch(`${BASE}/rest/v1/rpc/${fn}`, { method: 'POST', headers, body: JSON.stringify(args || {}) });
        const data = await res.json().catch(() => ({}));
        if (res.ok && data.ok) return data;
        lastErr = new Error(data.error || data.message || `Request failed (HTTP ${res.status})`);
        if (FATAL.test(lastErr.message) || (res.status >= 400 && res.status < 500 && res.status !== 429)) {
          lastErr.fatal = true;
          throw lastErr;
        }
      } catch (e) {
        lastErr = e;
        if (e.fatal) throw e;
      }
      if (i < tries - 1) await sleep(Math.min(4000, 500 * 2 ** i) * (0.5 + Math.random()));
    }
    throw lastErr;
  }

  const RPC = { questions: 'questions_for', results: 'results', moments: 'moments' };

  let memToken = null;
  function token() {
    if (memToken) return memToken;
    try { memToken = localStorage.getItem('lp_token'); } catch (e) { /* storage blocked */ }
    if (!memToken) {
      const raw = (window.crypto && crypto.randomUUID)
        ? crypto.randomUUID()
        : Math.random().toString(36).slice(2) + Date.now().toString(36);
      memToken = raw.replace(/[^A-Za-z0-9_-]/g, '');
      try { localStorage.setItem('lp_token', memToken); } catch (e) { /* fine */ }
    }
    return memToken;
  }

  function store(key, value) {
    try {
      if (value === undefined) return JSON.parse(localStorage.getItem(key) || 'null');
      localStorage.setItem(key, JSON.stringify(value));
    } catch (e) { return null; }
  }

  // Tiny DOM builder: h('div', {class: 'x', onclick: fn}, 'text', childNode)
  function h(tag, attrs, ...children) {
    const el = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs || {})) {
      if (v === null || v === undefined || v === false) continue;
      if (k.startsWith('on')) el.addEventListener(k.slice(2), v);
      else if (k === 'class') el.className = v;
      else if (k === 'style' && typeof v === 'object') Object.assign(el.style, v);
      else if (k in el && typeof v !== 'string') el[k] = v;
      else el.setAttribute(k, v === true ? '' : v);
    }
    for (const c of children.flat()) {
      if (c === null || c === undefined || c === false) continue;
      el.append(c instanceof Node ? c : document.createTextNode(String(c)));
    }
    return el;
  }

  window.LivePoll = {
    params, h, token, store,
    get: (action, query) => rpc(RPC[action], query && 'm' in query ? { p_moment: query.m } : {}),
    post: (body) => rpc('submit', { p_moment: body.moment, p_token: body.token, p_answers: body.answers }),
    configured: () => !!BASE && !BASE.includes('PASTE_') && (params.has('sb') || (!!KEY && !KEY.includes('PASTE_'))),
  };
})();
