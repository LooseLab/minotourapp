/**
 * Poll processing-queue API and update navbar badge (template_private).
 */
(function (window) {
  'use strict'

  function formatReads (n) {
    if (n == null || isNaN(n)) return '—'
    if (n >= 1000000) return (n / 1000000).toFixed(1).replace(/\.0$/, '') + 'M'
    if (n >= 10000) return Math.round(n / 1000) + 'k'
    if (n >= 1000) return (n / 1000).toFixed(1).replace(/\.0$/, '') + 'k'
    return String(n)
  }

  function render (el, data) {
    var textEl = el.querySelector('.obsidian-queue-indicator__text')
    if (!textEl) return

    if (data.error) {
      el.classList.remove(
        'obsidian-queue-indicator--loading',
        'obsidian-queue-indicator--idle',
        'obsidian-queue-indicator--busy'
      )
      el.classList.add('obsidian-queue-indicator--error')
      textEl.textContent = 'Queue ?'
      el.title = 'Could not load queue status (Redis unavailable?)'
      return
    }

    var batches = data.read_batches_queued
    var reads = data.reads_queued
    var harvesting = data.harvesting_active

    el.classList.remove(
      'obsidian-queue-indicator--loading',
      'obsidian-queue-indicator--error',
      'obsidian-queue-indicator--idle',
      'obsidian-queue-indicator--busy'
    )

    if (!batches && !reads) {
      el.classList.add('obsidian-queue-indicator--idle')
      textEl.textContent = 'Queue clear'
      el.title = 'No read batches waiting for harvest. Idle.'
      return
    }

    el.classList.add('obsidian-queue-indicator--busy')

    var parts = []
    if (reads != null && reads > 0) {
      parts.push(formatReads(reads) + ' reads')
    }
    if (batches != null && batches > 0) {
      parts.push(batches + (batches === 1 ? ' batch' : ' batches'))
    }
    textEl.textContent = parts.join(' · ') || 'Queued'

    var tip = []
    if (reads != null) tip.push(reads.toLocaleString() + ' reads')
    if (batches != null) tip.push(batches + ' batch(es) in Redis awaiting harvest')
    if (harvesting) tip.push('Harvest task running')
    el.title = tip.join('. ') + '.'
  }

  function fetchStatus (url, el) {
    if (typeof window.fetch !== 'function') return
    window.fetch(url, {
      credentials: 'same-origin',
      headers: { Accept: 'application/json' }
    }).then(function (r) {
      if (!r.ok) throw new Error('HTTP ' + r.status)
      return r.json()
    }).then(function (data) {
      render(el, data)
    }).catch(function () {
      render(el, { error: true })
    })
  }

  window.MinoTourQueueIndicator = {
    init: function (opts) {
      var el = document.getElementById('minotour-queue-indicator')
      if (!el) return
      var url = (opts && opts.url) || el.getAttribute('data-queue-url')
      if (!url) return
      var ms = (opts && opts.intervalMs) || 5000

      fetchStatus(url, el)
      window.setInterval(function () {
        fetchStatus(url, el)
      }, ms)
    }
  }
})(window)
