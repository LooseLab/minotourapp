const gulp = require('gulp');
const gutil = require('gulp-util');
const uglify = require('gulp-uglify');
const watch = require('gulp-watch');
const babel = require('gulp-babel');
const webpack = require('gulp-webpack');


gulp.task('scripts', function () {
    return gulp
        .src(['web/static/web/minotour.js'])
        .pipe(babel({
            presets: ['es2015']
        }))
        .pipe(gulp.dest('dist2'));
});

gulp.task('scripts2', function () {
    return gulp
        .src(['dist2/*.js'])
        .pipe(webpack())
        .pipe(gulp.dest('dist'));
});

gulp.task('watch', function () {
    var watcher = gulp.watch('web/static/web/minotour.js', ['scripts']);

    watcher.on('change', function (event) {
        gutil.log('File ' + event.path + ' was ' + event.type + ', running tasks...');
    });
});

gulp.task('teste', function () {
    return gulp
    // .src(['web/static/web/**/*.js'])
        .src(['web/static/web/minotour.js'])
        .pipe(webpack())
        .pipe(babel({
            presets: ['es2015']
        }))
        .pipe(uglify())
        .pipe(gulp.dest('dist'));
});