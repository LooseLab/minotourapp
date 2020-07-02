class MinotourTask {
  // a class for a Job to be created - stores all the attributes we might need
  constructor (flowcellId, jobTypeId, referenceId, taretSetId, fromDatabase) {
    // construct all the attributes we might need
    this._flowcellId = flowcellId
    this._jobTypeId = jobTypeId
    this._referenceId = referenceId
    this._targetSetId = taretSetId
    this._fromDatabase = fromDatabase
  }

  // getter for flowcell id
  get flowcellId () {
    return this._flowcellId
  }

  // getter for job type id
  get jobTypeId () {
    return this._jobTypeId
  }

  // getter for the reference ID
  get referenceId () {
    return this._referenceId
  }

  // getter for the target_set id
  get targetSetId () {
    return this._targetSetId
  }

  get fromDatabase () {
    return this._fromDatabase
  }
}
