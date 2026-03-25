/**
 * Settings drawer: AdminLTE ControlSidebar assumes .wrapper + .main-header; minoTour uses #tn only.
 * Inline styles on the aside also break slide open/close. We toggle body classes + optional backdrop only.
 */
(function ($) {
  'use strict'

  function isOpen () {
    return $('body').hasClass('control-sidebar-slide-open') || $('body').hasClass('control-sidebar-open')
  }

  function removeBackdrop () {
    $('#obsidian-settings-backdrop').remove()
  }

  function close () {
    $('body').removeClass('control-sidebar-slide-open control-sidebar-open')
    $('html').removeClass('control-sidebar-animate')
    removeBackdrop()
  }

  function openBackdrop () {
    if ($('#obsidian-settings-backdrop').length) return
    $('<div>', {
      id: 'obsidian-settings-backdrop',
      class: 'obsidian-settings-backdrop',
      'aria-hidden': 'true'
    }).appendTo('body').on('click', function () {
      close()
    })
  }

  function toggle (e) {
    if (e) {
      e.preventDefault()
      e.stopImmediatePropagation()
    }
    if (isOpen()) {
      close()
    } else {
      $('html').addClass('control-sidebar-animate')
      $('body').addClass('control-sidebar-slide-open')
      openBackdrop()
    }
  }

  $(function () {
    if (!$('.obsidian-settings-sidebar.control-sidebar').length) return

    $(document).off('click', '[data-widget="control-sidebar"]')
    $(document).on('click', '[data-widget="control-sidebar"]', toggle)

    $(document).on('click', '#minotour-settings-close', function (e) {
      e.preventDefault()
      close()
    })

    $(document).on('keydown.minotourSettings', function (e) {
      if (e.key === 'Escape' && isOpen()) {
        close()
      }
    })
  })
})(jQuery)
